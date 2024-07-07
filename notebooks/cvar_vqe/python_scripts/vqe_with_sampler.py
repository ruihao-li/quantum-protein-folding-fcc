#####

import sys

sys.path.append("../../../")


import json
import os
import time
from datetime import datetime, timezone

import execution_utils
import matplotlib.pyplot as plt
import numpy as np
import psutil
import ray
from execution_utils.measurements import (
    get_cvar_energy,
    process_counts_parallel,
    process_counts_serial,
)
from numpy.random import default_rng
from qiskit import QuantumCircuit
from qiskit.circuit.library import EfficientSU2, RealAmplitudes
from qiskit.quantum_info import SparsePauliOp
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_aer import AerSimulator
from qiskit_ibm_runtime import QiskitRuntimeService
from qiskit_ibm_runtime import SamplerV2 as Sampler
from qiskit_ibm_runtime import Session
from qufold import (
    MiyazawaJerniganInteraction,
    PenaltyParameters,
    Peptide,
    ProteinFoldingProblem,
)
from qufold.protein_folding_result import ProteinFoldingResult


NUM_WORKERS = (
    psutil.cpu_count()
)  # allowing Ray to use all core. set as per your choice/need
ray.init(
    num_cpus=NUM_WORKERS,
    log_to_driver=False,
    ignore_reinit_error=True,
    runtime_env={
        "py_modules": [execution_utils],
    },
)

TIMESTAMP = datetime.now(timezone.utc)


## build and set parameters for protein folding problem
def build_pf(main_seq: str):
    """Builds the protein folding problem for the given sequence."""
    # Define the interaction
    side_chains = [""] * len(main_seq)

    mj_interaction = MiyazawaJerniganInteraction()

    penalty_back = 100
    penalty_chiral = 100
    penalty_1 = 100

    penalty_terms = PenaltyParameters(penalty_chiral, penalty_back, penalty_1)

    peptide = Peptide(main_seq, side_chains)

    protein_folding_problem = ProteinFoldingProblem(
        peptide, mj_interaction, penalty_terms
    )

    return protein_folding_problem


"""Define problem parameters"""
# main chain
main_chain = "GSNQNNF"
# main_chain = "YYDPETGTWY"
# main_chain = "RGKWTYNGITYEGR"

# max number of bitstrings to save. Set `None` to save all.
MAX_NUM_SAVE_BITSTRING = 100_000

# check N lowest bitstrings.
# will save `.xyz` file for N lowest energy bitstrings in the corresponding directory
NUM_CHECK_LOWSEST_ENERGY_BITSRINGS = 50


"""Define protein folding problem"""
pf = build_pf(main_chain)  # creates the PF problem instance
qubit_op = pf.qubit_op()  # creates the problem Hamiltonian
print(f"Num qubits {qubit_op.num_qubits}")


"""Get backend"""
# hardware backend
# service = QiskitRuntimeService()
# backend = service.backend("ibm_cleveland")

# or simulator backend
backend = AerSimulator(method="matrix_product_state")

# or simulator from a hardware backend (noisy simulation)
# service = QiskitRuntimeService()
# backend = service.backend("ibm_cleveland")
# backend = AerSimulator.from_backend(backend=backend)


"""Define variational ansatz"""
# currently, uses either RealAmplitudes or EfficientSU2
# can be other ansatz, even custom ansatz of your choice
ansatz_type = "ra"
# ansatz_type = "su2"

if ansatz_type == "ra":
    ansatz = RealAmplitudes(
        num_qubits=qubit_op.num_qubits, reps=2, entanglement="pairwise"
    ).decompose() # `.decompose` is optional if you are transpiling
elif ansatz_type == "su2":
    ansatz = EfficientSU2(
        num_qubits=qubit_op.num_qubits, reps=2, entanglement="pairwise"
    ).decompose() # `.decompose` is optional if you are transpiling

# add measure ops for Sampler job
ansatz.measure_all()

# transpilation for HW runs. Optional for Simulator.
pm = generate_preset_pass_manager(backend=backend, optimization_level=3)
isa_circuit = pm.run(ansatz)


"""generate random initial parameter values for VQE loop"""
rng = default_rng(seed=0)
pi = np.pi
init_parameter_values = rng.uniform(-pi, pi, size=isa_circuit.num_parameters)


"""cost function"""
# global variables to save intermediate data
energies = []
parameters = []
bitstring_expval_all = {}
def cost_func(
    params: np.array,
    ansatz: QuantumCircuit,
    hamiltonian: SparsePauliOp,
    sampler: Sampler,
) -> float:
    """Return estimate of energy from estimator

    Args:
        params (ndarray): Array of ansatz parameters
        ansatz (QuantumCircuit): Parameterized ansatz circuit
        hamiltonian (SparsePauliOp): Operator representation of Hamiltonian
        sampler (SamplerV2): Estimator primitive instance

    Returns:
        float: Energy estimate
    """
    tic0 = time.time()
    pub = (ansatz, params)

    tic1 = time.time()
    primitive_result = sampler.run(pubs=[pub]).result()
    pub_result = primitive_result[0]
    counts = pub_result.data.meas.get_counts()
    toc1 = time.time()
    print(f" >> Sampler job took {toc1-tic1:.4f} seconds")

    tic2 = time.time()
    bitstring_wise_expval, prob_expval_list = process_counts_parallel(
        counts=counts, observable=hamiltonian, num_batches=NUM_WORKERS
    )
    bitstring_expval_all.update(bitstring_wise_expval)
    toc2 = time.time()
    print(
        f" >> processing {len(prob_expval_list)} unique bitstrings took {toc2-tic2:.4f} seconds"
    )

    energy = get_cvar_energy(measurements=prob_expval_list)
    print(f"energy {energy:.4f}\n")

    energies.append(energy)
    parameters.append(params)
    toc0 = time.time()
    print(f"One iteration of optimizer took {toc0-tic0:.4f} seconds")

    return energy


"""Execution using Runtime primitive"""
runner = "simulator" if backend.configuration().simulator else "hardware"
# Simulator runs does not require Session
# for HW runs, the optimization loop must be inside a Session
with Session(backend=backend) as session:
    sampler = Sampler(mode=session)
    sampler.options.default_shots = 10_000

    if runner == "hardware":
        # recommended: for Sampler runs enable gates twirling and disable measure twirling
        sampler.options.twirling.enable_gates = True
        sampler.options.twirling.enable_measure = False

        ## based on ansatz you may enable or disable dynamical decoupling (DD).
        ## for RealAmplitudes/EfficientSU2 with `pairwise` entanglement, disable DD
        # sampler.options.dynamical_decoupling.enable = True
        # sampler.options.dynamical_decoupling.sequence_type = "XX" # "XpXm" or "XY4"

    
    # use exactly one of the three optimizer blocks below.
    # NFT (from `qiskit_algorithms`) and COBYLA (from `scipy`) are example only
    # use any optimizer of your choice from `qiskit_algorithms` or `scipy`
    # set optimizer parameters such as maxiter as per your need
    """qiskit-algorithms optimizers"""
    # from qiskit_algorithms.optimizers import NFT, SPSA
    # optimizer = NFT(maxiter=600)
    # fun_wrap = NFT.wrap_function(cost_func, (isa_circuit, qubit_op, sampler))
    # optimizer_result = optimizer.minimize(fun_wrap, x0=init_parameter_values)
    # optimized_param_values = optimizer_result.x

    """CMA optimizer (https://github.com/CMA-ES/pycma) (pip install cma)"""
    import cma
    sigma0 = 0.5 # initial value for the variance for cma-es
    optimized_param_values, es = cma.fmin2(
        cost_func,
        init_parameter_values,
        sigma0,
        args=(isa_circuit, qubit_op, sampler),
        options={'maxiter': 100}
    )

    """SciPy optimizer"""
    # from scipy.optimize import minimize
    # optimizer_result = minimize(
    #     fun=cost_func,
    #     x0=init_parameter_values,
    #     method="cobyla",
    #     args=(isa_circuit, qubit_op, sampler),
    #     options={"maxiter": 5},
    # )
    # optimized_param_values = optimizer_result.x

print(optimized_param_values)


"""process and save results"""
timestamp_str = TIMESTAMP.strftime("%Y_%m_%d_%H_%M_%S_%f")
parent_dir = f"../data/{main_chain}_{timestamp_str}_{backend.name}"
xyz_files_dir = f"{parent_dir}/xyz_files"
plots_3d_dir = f"{parent_dir}/plots_3D_structure"

os.mkdir(parent_dir)
os.mkdir(xyz_files_dir)
os.mkdir(plots_3d_dir)

sorted_bitstring_expval_all = sorted(bitstring_expval_all.items(), key=lambda x: x[1])
with open(f"{parent_dir}/sorted_bitstrings_{main_chain}_{ansatz_type}.json", "w") as jf:
    json.dump(sorted_bitstring_expval_all[:MAX_NUM_SAVE_BITSTRING], jf, indent=2)

# plot optimization convergance
fig = plt.figure()

plt.plot(range(len(energies)), np.real(energies))
plt.ylabel("Conformational Energy")
plt.xlabel("VQE Iterations")

# fig.add_axes([0.44, 0.51, 0.44, 0.32])
# plt.plot(counts[40:], values[40:])
# plt.ylabel("Conformation Energy")
# plt.xlabel("VQE Iterations")

plt.savefig(
    f"{parent_dir}/conf_energy_plot_{runner}_{main_chain}.png",
    dpi=300,
    transparent=False,
)

for idx in range(NUM_CHECK_LOWSEST_ENERGY_BITSRINGS):
    # manually creating `ProteinFoldingResult` with the bitstring with lowest energies
    # sorted_bitstring_expval_all: List[Tuple[str, float]]
    # 1st elem is the bitstring, 2nd elem is energy
    result = ProteinFoldingResult(
        unused_qubits=pf.unused_qubits,
        peptide=pf.peptide,
        turn_sequence=sorted_bitstring_expval_all[idx][0],
    )

    # optinal lines to print extra information
    # print(
    #     "The bitstring representing the shape of the protein during optimization is: ",
    #     result.turn_sequence,
    # )
    # print("The expanded expression is:", result.get_result_binary_vector())

    # print(f"The folded protein's main sequence of turns is: {result.protein_shape_decoder.main_turns}")
    # print(f"and the side turn sequences are: {result.protein_shape_decoder.side_turns}")

    # save 3D structures of the protein
    fig = result.get_figure(title="3dcrd", ticks=False, grid=True)
    fig.get_axes()[0].view_init(10, 70)
    fig.savefig(
        f"{plots_3d_dir}/{idx:04}_3Dcrd_{runner}_{main_chain}.png",
        dpi=300,
        transparent=False,
    )
    # plt.show()
    plt.close()

    # saves xyz files
    result.save_xyz_file(
        name=f"{idx:04}_{runner}_{main_chain}",
        path=xyz_files_dir,
        comment=f"Generated using {backend.name}. Bitstring '{sorted_bitstring_expval_all[idx][0]}'",
        replace=True,
    )
