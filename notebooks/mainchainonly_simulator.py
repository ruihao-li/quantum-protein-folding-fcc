## this is part of the original code, using Qiskit Nature. For simulator jobs, you can use the
## version now available on Qiskit Research. For hardware jobs however, we are still trying to figure out
## how to submit when the code is using Qiskit Research.

from qiskit_nature.problems.sampling.protein_folding.interactions.random_interaction import (
    RandomInteraction,
)
from qiskit_nature.problems.sampling.protein_folding.interactions.miyazawa_jernigan_interaction import (
    MiyazawaJerniganInteraction,
)
from qiskit_nature.problems.sampling.protein_folding.peptide.peptide import Peptide
from qiskit_nature.problems.sampling.protein_folding.protein_folding_problem import (
    ProteinFoldingProblem,
)

from qiskit_nature.problems.sampling.protein_folding.penalty_parameters import PenaltyParameters

from qiskit.utils import algorithm_globals, QuantumInstance

algorithm_globals.random_seed = 23

##

main_chain = "LHPGAGK"
side_chains = [""] * len(main_chain)

##

random_interaction = RandomInteraction()
mj_interaction = MiyazawaJerniganInteraction()

##

penalty_back = 10
penalty_chiral = 10
penalty_1 = 10

penalty_terms = PenaltyParameters(penalty_chiral, penalty_back, penalty_1)
##
peptide = Peptide(main_chain, side_chains)

##

protein_folding_problem = ProteinFoldingProblem(peptide, mj_interaction, penalty_terms)
qubit_op = protein_folding_problem.qubit_op()

with open("hamiltonian.txt", "w+") as hamiltonian_file:
     hamiltonian_file.write(str(qubit_op))

##

print(qubit_op.num_qubits)

##

from qiskit.circuit.library import RealAmplitudes
from qiskit.algorithms.optimizers import COBYLA
from qiskit.algorithms import NumPyMinimumEigensolver, VQE
from qiskit.opflow import PauliExpectation, CVaRExpectation
from qiskit import execute, Aer
from qiskit.providers.ibmq import least_busy
from qiskit import IBMQ
IBMQ.load_account()
IBMQ.providers()
provider = IBMQ.get_provider(hub='ibm-q-ccf', group='discovery-accel', project='qperf')
small_devices = provider.backends(filters=lambda x: x.configuration().n_qubits >= qubit_op.num_qubits
                                   and not x.configuration().simulator)
# set classical optimizer
optimizer = COBYLA(maxiter=50)

# set variational ansatz
ansatz = RealAmplitudes(reps=1)

# set the backend
backend_name = "aer_simulator"
# backend = least_busy(small_devices)
backend = QuantumInstance(
    Aer.get_backend(backend_name),
    shots=8192,
    seed_transpiler=algorithm_globals.random_seed,
    seed_simulator=algorithm_globals.random_seed,
)
print(backend)

counts = []
values = []

def store_intermediate_result(eval_count, parameters, mean, std):
    counts.append(eval_count)
    values.append(mean)


# initialize CVaR_alpha objective with alpha = 0.1
cvar_exp = CVaRExpectation(0.1, PauliExpectation())

# initialize VQE using CVaR
vqe = VQE(
    expectation=cvar_exp,
    optimizer=optimizer,
    ansatz=ansatz,
    quantum_instance=backend,
    callback=store_intermediate_result,
)

raw_result = vqe.compute_minimum_eigenvalue(qubit_op)
print(raw_result)
with open("raw_result.txt", "w+") as rawresults_file:
     rawresults_file.write(str(raw_result))
##

with open("counts.txt", "w+") as countsfile, open("values.txt", "w+") as valuesfile:
     for num in counts:
         countsfile.write("{}\n".format(num))
     for num in values:
         valuesfile.write("{}\n".format(num))

import matplotlib.pyplot as plt

fig = plt.figure()

plt.plot(counts, values)
plt.ylabel("Conformation Energy")
plt.xlabel("VQE Iterations")

fig.add_axes([0.44, 0.51, 0.44, 0.32])

plt.plot(counts[40:], values[40:])
plt.ylabel("Conformation Energy")
plt.xlabel("VQE Iterations")
plt.savefig("conf_energy_simulator.png", dpi=300,transparent=False)
#plt.show()
##

result = protein_folding_problem.interpret(raw_result=raw_result)
print(
    "The bitstring representing the shape of the protein during optimization is: ",
    result.turn_sequence,
)
print("The expanded expression is:", result.get_result_binary_vector())

##

print(f"The folded protein's main sequence of turns is: {result.protein_shape_decoder.main_turns}")
print(f"and the side turn sequences are: {result.protein_shape_decoder.side_turns}")


result.save_xyz_file(replace=True)

##

fig = result.get_figure(title="3dcrd", ticks=False, grid=True)
fig.get_axes()[0].view_init(10, 70)
#fig.show()
fig.savefig("3dcrd_Simulator.png", dpi=300,transparent=False)
##


import qiskit.tools.jupyter


#%qiskit_version_table
#%qiskit_copyright
