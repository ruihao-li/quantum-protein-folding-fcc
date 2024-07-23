#####

import sys

sys.path.append("../../../")


import json
import os
import time
from datetime import datetime, timezone

import matplotlib.pyplot as plt
import numpy as np
import psutil
import qufold
import ray
from numpy.random import default_rng
from qiskit import QuantumCircuit
from qiskit.circuit.library import EfficientSU2, RealAmplitudes
from qiskit.quantum_info import SparsePauliOp
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_aer import AerSimulator
from qiskit_ibm_runtime import QiskitRuntimeService, RuntimeEncoder
from qiskit_ibm_runtime import SamplerV2 as Sampler
from qiskit_ibm_runtime import Session
from qufold import (
    MiyazawaJerniganInteraction,
    PenaltyParameters,
    Peptide,
    ProteinFoldingProblem,
)
from qufold.execution_utils.measurements import (
    get_cvar_energy,
    process_counts_parallel,
    process_counts_serial,
)
from qufold.protein_folding_result import ProteinFoldingResult

NUM_WORKERS = 10  # allowing Ray to use all core. set as per your choice/need
ray.init(
    num_cpus=NUM_WORKERS,
    log_to_driver=False,
    ignore_reinit_error=True,
    runtime_env={
        "py_modules": [qufold],
    },
)

metadata = {}
TIMESTAMP = datetime.now(timezone.utc)

metadata["timestamp"] = TIMESTAMP.isoformat()

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
# main_chain = "YYDPETGTWY" #PDB: 5AWL chignolin, 10 AA
# main_chain = "FATMRYPSDSDE" #PDB: 1IXU, 12 AA
main_chain = "VRRFDLLKRILK" #PDB: 2N5R, 12 AA
# main_chain = "RGKWTYNGITYEGR" #PDB:1K43
# main_chain = "VLAMWKVGFFKRNRP" # Jun's 15 AA structure
# main_chain = "GGLRSLGRKILRAWKKYG" #PDB:2NDC, 18AA
# main_chain = "IGLRGLGRKIALIHKKYG" #PDB:2NDE, 18AA
# main_chain = "DAYAQWLKDGGPSSGRPPPS" #trp-cage, 20 AA
# main_chain = "GNDYEDRYYRENMYRYPNQVYYRPVC" #PDB:1G04, 26AA
# main_chain = "YYHFWHRGVTKRSLSPHRPRHSRLQR" #PDB:6A8Y, 26AA

# shots for circuit execution
SHOTS = 100_000

# max number of bitstrings to save. Set `None` to save all.
MAX_NUM_SAVE_BITSTRING = 100_000

# check N lowest bitstrings.
# will save `.xyz` file for N lowest energy bitstrings in the corresponding directory
NUM_CHECK_LOWSEST_ENERGY_BITSRINGS = 20

metadata["shots"] = SHOTS
metadata["num_workers"] = NUM_WORKERS
metadata["main_chain"] = main_chain
metadata["max_num_save_bitstring"] = MAX_NUM_SAVE_BITSTRING
metadata["num_check_lowest_energy_bitstrings"] = NUM_CHECK_LOWSEST_ENERGY_BITSRINGS

"""Define protein folding problem"""
pf = build_pf(main_chain)  # creates the PF problem instance
qubit_op: SparsePauliOp = pf.qubit_op()  # creates the problem Hamiltonian
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

metadata["backend_name"] = backend.name

"""Define variational ansatz"""
# currently, uses either RealAmplitudes or EfficientSU2
# can be other ansatz, even custom ansatz of your choice
ansatz_type = "ra"
# ansatz_type = "su2"

if ansatz_type == "ra":
    ansatz = RealAmplitudes(
        num_qubits=qubit_op.num_qubits, reps=2, entanglement="pairwise"
    ).decompose()  # `.decompose` is optional if you are transpiling
elif ansatz_type == "su2":
    ansatz = EfficientSU2(
        num_qubits=qubit_op.num_qubits, reps=2, entanglement="pairwise"
    ).decompose()  # `.decompose` is optional if you are transpiling

# add measure ops for Sampler job
ansatz.measure_all()

# transpilation for HW runs. Optional for Simulator.
pm = generate_preset_pass_manager(backend=backend, optimization_level=3)
isa_circuit = pm.run(ansatz)

metadata["physical_circuit"] = RuntimeEncoder().encode(isa_circuit)

"""generate random initial parameter values for VQE loop"""
rng = default_rng(seed=0)
pi = np.pi
init_parameter_values = rng.uniform(-pi, pi, size=isa_circuit.num_parameters)

metadata["initial_param_values"] = init_parameter_values.tolist()

"""cost function"""
# global variables to save intermediate data
energies = []
parameters = []
job_ids = []
func_eval_durations = []
primitive_job_durations = []
bs_processing_durations = []
bitstring_expval_all = {}
nfev = 0

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
    global nfev
    
    tic0 = time.time()
    pub = (ansatz, params)

    tic1 = datetime.now(timezone.utc)
    print(f"# cost function eval: {nfev+1:04}")
    print(f"Submitting Sampler job at time {tic1.isoformat()}")
    # submitting sampler job in try-except block.
    # If there is an error submitting the job (e.g., Session is closed)
    # Exception block will capture it, and return 'inf'
    try:
        job = sampler.run(pubs=[pub])
        job_id = job.job_id()
        print(f"Job ID: {job_id}")
        primitive_result = job.result()
    except Exception as e:
        print("Exception happened during Runtime Sampler job")
        print(f"Session status {sampler.session.status()}")
        print(repr(e))
        print("\n")
        return float("inf")
    pub_result = primitive_result[0]
    counts = pub_result.data.meas.get_counts()
    toc1 = datetime.now(timezone.utc)
    job_duration = (toc1 - tic1).total_seconds()

    job_ids.append(job_id)
    
    tic2 = time.time()
    bitstring_wise_expval, prob_expval_list = process_counts_parallel(
        counts=counts, observable=hamiltonian, num_batches=NUM_WORKERS
    )
    bitstring_expval_all.update(bitstring_wise_expval)
    energy = get_cvar_energy(measurements=prob_expval_list)
    toc2 = time.time()
    
    energies.append(energy)
    parameters.append(params)
    toc0 = time.time()
    
    fev_dur = toc0-tic0
    bs_processing_dur = toc2-tic2
    func_eval_durations.append(round(fev_dur, 4))
    bs_processing_durations.append(round(bs_processing_dur, 4))
    primitive_job_durations.append(round(job_duration, 4))
    
    print(f"One iteration of cost func took {fev_dur:.4f} seconds")
    print(f" >> Sampler job took {job_duration:.4f} seconds")
    print(
        f" >> processing {len(prob_expval_list)}"
        f" unique bitstrings took {bs_processing_dur:.4f} seconds"
    )
    print(f" >> energy {energy:.4f}\n")
    nfev += 1
    return energy


"""Execution using Runtime primitive"""
runner = "simulator" if backend.configuration().simulator else "hardware"
# Simulator runs does not require Session
# for HW runs, the optimization loop must be inside a Session
with Session(backend=backend) as session:
    metadata["session_id"] = session.session_id
    sampler = Sampler(mode=session)
    sampler.options.default_shots = SHOTS

    if runner == "hardware":
        # recommended: for Sampler runs enable gates twirling and disable measure twirling
        sampler.options.twirling.enable_gates = False
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
    # import cma

    # sigma0 = 0.5  # initial value for the variance for cma-es
    # optimized_param_values, es = cma.fmin2(
    #     cost_func,
    #     init_parameter_values,
    #     sigma0,
    #     args=(isa_circuit, qubit_op, sampler),
    #     options={"maxiter": 3},
    # )

    """SciPy optimizer"""
    from scipy.optimize import minimize
    optimizer_result = minimize(
        fun=cost_func,
        x0=init_parameter_values,
        method="cobyla",
        args=(isa_circuit, qubit_op, sampler),
        options={"maxiter": 10},
    )
    optimized_param_values = optimizer_result.x

print(f"optimized parameter values: {optimized_param_values}")

metadata["optimized_param_values"] = optimized_param_values.tolist()
metadata["nfev"] = nfev
metadata["job_ids"] = job_ids
metadata["durations (sec)"] = {
    "func_eval_durations": func_eval_durations,
    "job_durations": primitive_job_durations,
    "bitstring_processing_durations": bs_processing_durations,
}

"""process and save results"""
timestamp_str = TIMESTAMP.strftime("%Y_%m_%d_%H_%M_%S_%f")
parent_dir = f"../data/{main_chain}_{timestamp_str}_{backend.name}"
xyz_files_dir = f"{parent_dir}/xyz_files"
plots_3d_dir = f"{parent_dir}/plots_3D_structure"

os.mkdir(parent_dir)
os.mkdir(xyz_files_dir)
os.mkdir(plots_3d_dir)

with open(f"{parent_dir}/metadata.json", "w") as jf:
    json.dump(metadata, jf, indent=2, sort_keys=True)

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
