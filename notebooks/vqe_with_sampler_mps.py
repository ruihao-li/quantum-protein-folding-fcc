#####

import sys

sys.path.append("../../../../")
sys.path.append("../")

from qiskit.circuit.library import RealAmplitudes, EfficientSU2
from qiskit_ibm_runtime import SamplerV2 as Sampler

import matplotlib.pyplot as plt
import numpy as np
from qufold import (
    MiyazawaJerniganInteraction,
    Peptide,
    ProteinFoldingProblem,
    PenaltyParameters,
)
from qiskit.quantum_info import SparsePauliOp

#####
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

    protein_folding_problem = ProteinFoldingProblem(peptide, mj_interaction, penalty_terms)

    return protein_folding_problem

#####

main_chain = "GSNQNNF"
# main_chain = "YYDPETGTWY"
# main_chain = "RGKWTYNGITYEGR"
pf = build_pf(main_chain) #creates the PF problem instance

#####

qubit_op = pf.qubit_op() #creates the problem Hamiltonian
print(f"Num qubits {qubit_op.num_qubits}")
#####

from qiskit_ibm_runtime import Session, Options, QiskitRuntimeService, SamplerV2 as Sampler
from qiskit_aer import AerSimulator

# service = QiskitRuntimeService()
# backend = service.backend("ibm_cleveland")
backend = AerSimulator(method="matrix_product_state")

#####
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

# set variational ansatz
# ansatz_type = "ra"
ansatz_type = "su2"
if ansatz_type == "ra":
    ansatz = RealAmplitudes(num_qubits=qubit_op.num_qubits, reps=3, entanglement="pairwise").decompose()
elif ansatz_type == "su2":
    ansatz = EfficientSU2(num_qubits=qubit_op.num_qubits, reps=2, entanglement="pairwise").decompose()
else:
    raise ValueError(f"'ansatz_type' {ansatz_type} must be 'ra' or 'su2'")

# add measure ops for Sampler job
ansatz_with_meas = ansatz.measure_all(inplace=False)


#####
import json
import numpy as np
from numpy.random import default_rng

rng = default_rng(seed=0)
pi = np.pi

init_parameter_values = rng.uniform(-pi, pi, size=ansatz_with_meas.num_parameters)

#####
import time
from cvar import get_cvar_energy, process_counts, process_counts_serial

counts_all = []
energies = []
parameters = []
bitstring_expval_all = {}

def cost_func(params, ansatz, hamiltonian, sampler: Sampler):
    """Return estimate of energy from estimator

    Parameters:
        params (ndarray): Array of ansatz parameters
        ansatz (QuantumCircuit): Parameterized ansatz circuit
        hamiltonian (SparsePauliOp): Operator representation of Hamiltonian
        sampler (SamplerV2): Estimator primitive instance

    Returns:
        float: Energy estimate
    """
    tic = time.time()
    pub = (ansatz, params)
    
    tic2 = time.time()
    primitive_result = sampler.run(pubs=[pub]).result()
    pub_result = primitive_result[0]
    counts = pub_result.data.meas.get_counts()
    toc2 = time.time()
    print(f' >> Simulation took {toc2-tic2:.4f} seconds')
    # counts_all.append(counts)
    
    tic1 = time.time()
    bitstring_wise_expval, prob_expval_list = process_counts(counts=counts, observable=hamiltonian)
    bitstring_expval_all.update(bitstring_wise_expval)
    toc1 = time.time()
    print(f' >> processing {len(prob_expval_list)} unique bitstrings took {toc1-tic1:.4f} seconds')

    energy = get_cvar_energy(measurements=prob_expval_list)
    
    print(f'energy {energy:.4f}\n')
    
    energies.append(energy)
    # parameters.append(params)
    toc = time.time()
    print(f'One Sampler call took {toc-tic:.4f} seconds')

    return energy

#####
from scipy.optimize import minimize
# import cma

# Simulator runs does not require Session
# for HW runs, move the optimization loop inside a Session
# also simulator runs do not have notion of error suppression and mitigation
# (e.g., dynamical decoupling, twirling, etc.). No need to set those options
sampler = Sampler(mode=backend)
sampler.options.default_shots = 10_000


"""qiskit-algorithms optimizers"""
# from qiskit_algorithms.optimizers import COBYLA, NFT

# optimizer = NFT(maxiter=600)
# fun_wrap = NFT.wrap_function(cost_func, (ansatz_with_meas, qubit_op, sampler))
# optimizer_result = optimizer.minimize(fun_wrap, x0=init_parameter_values)
# optimized_param_values = optimizer_result.x


"""CMA optimizer (https://github.com/CMA-ES/pycma)"""
import cma
sigma0 = 0.5 # initial value for the variance for cma-es
optimized_param_values, es = cma.fmin2(
    cost_func,
    init_parameter_values,
    sigma0, 
    args=(ansatz_with_meas, qubit_op, sampler),
    options={'maxiter': 10}
)
print(optimized_param_values)

sorted_bitstring_expval_all = sorted(bitstring_expval_all.items(), key=lambda x:x[1])
with open(f"sorted_bitstrings_{main_chain}_{ansatz_type}.json", "w") as jf:
    json.dump(sorted_bitstring_expval_all, jf, indent=2)


#####

import matplotlib.pyplot as plt

fig = plt.figure()

plt.plot(range(len(energies)), np.real(energies))
plt.ylabel("Conformational Energy")
plt.xlabel("VQE Iterations")

# fig.add_axes([0.44, 0.51, 0.44, 0.32])

# plt.plot(counts[40:], values[40:])
# plt.ylabel("Conformation Energy")
# plt.xlabel("VQE Iterations")
plt.savefig("conf_energy_plot_hardware_{}.png".format(main_chain), dpi=300,transparent=False)
plt.show()
# plt.close()
#####

from qufold.protein_folding_result import ProteinFoldingResult

# manually creating `ProteinFoldingResult` with the bitstring with lowest energy
# it is the first bitstring in `sorted_bitstring_expval_all`
# ideally you should check many candidate bitstrings instead of just one
# this snippet is for example only
# sorted_bitstring_expval_all: List[Tuple[str, float]] - 1st elem is the bitstring, 2nd elem is energy
result = ProteinFoldingResult(
    unused_qubits=pf.unused_qubits,
    peptide=pf.peptide,
    turn_sequence=sorted_bitstring_expval_all[0][0],
)

print(
    "The bitstring representing the shape of the protein during optimization is: ",
    result.turn_sequence,
)
print("The expanded expression is:", result.get_result_binary_vector())


##

print(f"The folded protein's main sequence of turns is: {result.protein_shape_decoder.main_turns}")
print(f"and the side turn sequences are: {result.protein_shape_decoder.side_turns}")

fig = result.get_figure(title="3dcrd", ticks=False, grid=True)
fig.get_axes()[0].view_init(10, 70)
fig.savefig("3Dcrd_hardware_{}.png".format(main_chain), dpi=300,transparent=False)
plt.show()
#####

result.save_xyz_file(replace=True)
