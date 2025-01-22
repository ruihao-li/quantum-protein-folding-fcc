iimport sys

sys.path.append("../")
from fcc import (
    MiyazawaJerniganInteraction,
    Peptide,
    ProteinFoldingProblem,
    PenaltyParameters,
    ProteinSolver,
)
import fcc
from qiskit.circuit.library import RealAmplitudes

from qiskit_ibm_runtime import SamplerV2 as Sampler
from qiskit_ibm_runtime.fake_provider import FakeSherbrooke  # for local testing mode
from qiskit_ibm_runtime import QiskitRuntimeService, Session
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
import ray
import psutil
import json
from time import time
import numpy as np
from datetime import datetime, timezone
import os


metadata = {}

# Define parameters
MAIN_SEQ = "GNLVS"
NUM_WORKERS = None
R2_THRESHOLD = 0.99
SHOTS = 10_000
OPTIMIZER = "COBYLA"
MAX_ITER = 200
MAX_NUM_SAVED_STATES = 1000
TIMESTAMP = datetime.now(timezone.utc)

metadata["main_seq"] = MAIN_SEQ
metadata["r2_threshold"] = R2_THRESHOLD
metadata["shots"] = SHOTS
metadata["optimizer"] = OPTIMIZER
metadata["max_iter"] = MAX_ITER
metadata["max_num_saved_states"] = MAX_NUM_SAVED_STATES
metadata["timestamp"] = TIMESTAMP.isoformat()

####################
num_workers = (
    NUM_WORKERS if NUM_WORKERS is not None else psutil.cpu_count(logical=False)
)
print(f"Number of workers: {num_workers}")
metadata["num_workers"] = num_workers
ray.init(
    num_cpus=num_workers,
    ignore_reinit_error=True,
    log_to_driver=False,
    runtime_env={
        "py_modules": [fcc],
    },
)


def build_pf(main_seq: str, energy_matrix_file: str = "mj_matrix"):
    """Builds the protein folding problem for the given sequence."""

    mj_interaction = MiyazawaJerniganInteraction(energy_matrix_file)
    # print(mj_interaction.calculate_energy_matrix(main_seq))

    penalty_back = 50
    penalty_redun = 50
    penalty_olap = 50

    penalty_terms = PenaltyParameters(penalty_back, penalty_redun, penalty_olap)

    peptide = Peptide(main_seq)

    protein_folding_problem = ProteinFoldingProblem(
        peptide, mj_interaction, penalty_terms
    )

    return protein_folding_problem


pf_problem = build_pf(MAIN_SEQ)
time_start = time()
qubit_op = pf_problem.qubit_op(r2_threshold=R2_THRESHOLD)
time_end = time()
print(f"Number of qubits: {qubit_op.num_qubits}")
print(f"Number of Hamiltonian terms: {len(qubit_op)}")
print(f"Time to generate the Hamiltonian operator: {time_end - time_start:.2f}s")

metadata["num_qubits"] = qubit_op.num_qubits
metadata["num_ham_terms"] = len(qubit_op)
metadata["ham_gen_time (s)"] = np.round(time_end - time_start, 2)


ansatz = RealAmplitudes(qubit_op.num_qubits, reps=1).decompose()
# Add measurements
ansatz.measure_all()
# Get backend
time_start = time()
service = QiskitRuntimeService()
backend = service.backend("ibm_cleveland")
# backend = FakeSherbrooke()
print(f"QPU backend: {backend.name}")
metadata["backend"] = backend.name
# Transpilation for hardware runs
pass_manager = generate_preset_pass_manager(backend=backend, optimization_level=3)
isa_circ = pass_manager.run(ansatz)
time_end = time()
print(f"Time to prepare the circuit: {time_end - time_start:.2f}s")
metadata["circ_preparation_time (s)"] = np.round(time_end - time_start, 2)

with Session(backend=backend) as session:
    metadata["session_id"] = session.session_id
    sampler = Sampler(mode=session)
    sampler.options.default_shots = SHOTS
    # recommended: for Sampler runs enable gates twirling and disable measure twirling
    sampler.options.twirling.enable_gates = False
    sampler.options.twirling.enable_measure = False
    protein_solver = ProteinSolver(
        ansatz=isa_circ, hamiltonian=qubit_op, sampler=sampler
    )
    opt_results = protein_solver.train(
        optimizer=OPTIMIZER,
        maxiter=MAX_ITER,
        num_batches=num_workers,
        num_saved_states=MAX_NUM_SAVED_STATES,
        verbose=True,
    )
print("VQE finished.")

# Save results and metadata
timestamp_str = TIMESTAMP.strftime("%Y-%m-%d")
parent_dir = f"../res/fcc_hw/{MAIN_SEQ}_{timestamp_str}_{backend.name}"
try:
    os.makedirs(parent_dir)
except FileExistsError:
    pass

with open(f"{parent_dir}/metadata.json", "w") as jf:
    json.dump(metadata, jf, indent=2, sort_keys=True)

with open(f"{parent_dir}/opt_results.json", "w") as jf:
    json.dump(opt_results, jf, indent=2)

print("Results saved.")

ray.shutdown()
