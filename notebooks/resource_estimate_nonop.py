## iteratively estimates resource requirement per sequence length, without optimized circuits

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

from qiskit.circuit.library import RealAmplitudes

from qiskit.compiler import transpile

from time import time

algorithm_globals.random_seed = 23

from qiskit.circuit.library import RealAmplitudes
from qiskit.algorithms.optimizers import COBYLA
from qiskit.algorithms import NumPyMinimumEigensolver, VQE
from qiskit.opflow import PauliExpectation, CVaRExpectation
from qiskit.providers.ibmq import least_busy
from qiskit import execute, Aer
from qiskit import IBMQ, Aer

provider = IBMQ.load_account()

IBMQ.providers()

provider = IBMQ.get_provider(hub='ibm-q-ccf', group='discovery-accel', project='qperf')

provider.backends()

# main_chain_big = "VLAMWKVGFFKRNRPVLAMWKV"
main_chain_big = "VLAMWKV"

with open("estimate_nonop.tsv", "w+") as outfile:
    outfile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format("size (AAs)", "qubits", "gate types", "depth", "optimized ansatz level", "gate type"))
    t0 = time()
    for i in range(4,len(main_chain_big)):
        print("index", i)
        main_chain = main_chain_big[0:i+1]
        print("main_chain", main_chain)
        side_chains = [""] * len(main_chain)

        random_interaction = RandomInteraction()
        mj_interaction = MiyazawaJerniganInteraction()

        penalty_back = 10
        penalty_chiral = 10
        penalty_1 = 10
        penalty_terms = PenaltyParameters(penalty_chiral, penalty_back, penalty_1)

        peptide = Peptide(main_chain, side_chains)

        protein_folding_problem = ProteinFoldingProblem(peptide, mj_interaction, penalty_terms)
        qubit_op = protein_folding_problem.qubit_op()
        qubit_op.num_qubits
        print("Qubit #", qubit_op.num_qubits)

        #define small devices based on qubit_op
#         small_devices = provider.backends(filters=lambda x: x.configuration().n_qubits >= qubit_op.num_qubits and not x.configuration().simulator)

        instructions = qubit_op.exp_i()
        print(type(instructions))

        # set classical optimizer
        optimizer = COBYLA(maxiter=50)

        # set variational ansatz
        # ansatz = RealAmplitudes(reps=1)
#         ansatz = RealAmplitudes(num_qubits=qubit_op.num_qubits, reps=3)
#         print(ansatz)
        # set the backend
        backend_name = "ibm_washington"
        # backend = least_busy(small_devices)
        backend = QuantumInstance(
            provider.get_backend(backend_name),
            shots=8192,
            seed_transpiler=algorithm_globals.random_seed,
            seed_simulator=algorithm_globals.random_seed,
        )

        counts = []
        values = []


        def store_intermediate_result(eval_count, parameters, mean, std):
            counts.append(eval_count)
            values.append(mean)


        ansatz = RealAmplitudes(num_qubits=qubit_op.num_qubits, reps=3)
#         print("ansatz")
#         print(ansatz)
        ansatz.count_ops()
        print("ansatz ops", ansatz.count_ops())
        print("gates", ansatz.count_ops().keys())
        print("gates #", len(ansatz.count_ops()))
        print("depth", ansatz.depth())
        
        optimized_0 = transpile(ansatz, backend=provider.get_backend(backend_name), seed_transpiler=11, optimization_level=0)            
        
        #iteratively write the resources required per sequence size 
        outfile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format((len(main_chain)), (qubit_op.num_qubits), (len(optimized_0.count_ops())), (optimized_0.depth()), "0", optimized_0.count_ops()))
    t1 = time()
    print("Total elapsed time {}s".format(int(t1 - t0)))


import qiskit.tools.jupyter

# qiskit_version_table
# qiskit_copyright
