import sys
import os
import json
import ray
import psutil
import numpy as np
from math import pi
from datetime import datetime, timezone

# Add repo root to path
sys.path.append("../../")

from fcc import (
    MiyazawaJerniganInteraction,
    Peptide,
    ProteinFoldingProblem,
    PenaltyParameters,
    ProteinSolver,
)
import fcc

from qiskit.circuit.library import RealAmplitudes
from qiskit.quantum_info import SparsePauliOp
from qiskit_ibm_runtime import QiskitRuntimeService, Session, SamplerV2 as Sampler
from qiskit_aer import AerSimulator
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

# === Metadata & Parameters ===
metadata = {}
MAIN_SEQ = "KLVFFA"
NUM_WORKERS = None
R2_THRESHOLD = 1.0
CHUNK = None
ANSATZ_REPS = 2
RUNNER = "aer-mps"
PARALLELIZER = "ray"
SHOTS = 100_000
OPTIMIZER = "COBYLA"
MAX_ITER = 500
MAX_NUM_SAVED_STATES = 1000
TIMESTAMP = datetime.now(timezone.utc)
metadata.update({
    "main_seq": MAIN_SEQ,
    "r2_threshold": R2_THRESHOLD,
    "chunk_size": CHUNK,
    "ansatz_reps": ANSATZ_REPS,
    "runner": RUNNER,
    "parallelizer": PARALLELIZER,
    "shots": SHOTS,
    "optimizer": OPTIMIZER,
    "max_iter": MAX_ITER,
    "max_num_saved_states": MAX_NUM_SAVED_STATES,
    "timestamp": TIMESTAMP.isoformat(),
})

# === Load Hamiltonian ===
# Use the full absolute path directly — don't redefine ham_path using os.path.join
ham_path = "/home/hakan-doga/data/protein-folding-qc/fcc_hamilts/KLVFFA_6AAs_R2_1000.json"

with open(ham_path, "r") as f:
    data = json.load(f)
    num_qubits = data["num_qubits"]
    sparse_op_list = [(label, coeff) for label, coeff in data["sparse_list"]]

qubit_op = SparsePauliOp.from_list(sparse_op_list, num_qubits=num_qubits)
print("qubit_op loaded")

metadata["num_qubits"] = qubit_op.num_qubits
metadata["num_ham_terms"] = len(qubit_op)

# === Prepare ansatz ===
ansatz = RealAmplitudes(qubit_op.num_qubits, reps=ANSATZ_REPS).decompose()
ansatz.measure_all()

if RUNNER == "aer-sv":
    backend = AerSimulator(method="statevector")
elif RUNNER == "aer-mps":
    backend = AerSimulator(method="matrix_product_state")
elif RUNNER == "hardware":
    service = QiskitRuntimeService()
    backend = service.backend("ibm_cleveland")

metadata["backend"] = backend.name
pass_manager = generate_preset_pass_manager(backend=backend, optimization_level=3)
isa_circ = pass_manager.run(ansatz)

# === Parallel Setup ===
num_workers = NUM_WORKERS or psutil.cpu_count(logical=False)
metadata["num_workers"] = num_workers
if PARALLELIZER == "ray":
    ray.init(
        num_cpus=num_workers,
        ignore_reinit_error=True,
        log_to_driver=False,
        runtime_env={"py_modules": [fcc]},
    )

# === Save Directory ===
base_save_dir = "/home/hakan-doga/data/protein-folding-qc/res/cheb_mps"
os.makedirs(base_save_dir, exist_ok=True)


# === VQE Loop ===
for run_id in range(10):
    print(f"\n--- VQE Run {run_id+1}/10 with Init Range: [4.71, 6.28] ---")
    init_params = np.random.uniform(3 * pi / 2, 2 * pi, size=isa_circ.num_parameters)

    with Session(backend=backend) as session:
            metadata["session_id"] = session.session_id
            sampler = Sampler(mode=session)
            sampler.options.default_shots = SHOTS
            if RUNNER == "hardware":
                sampler.options.twirling.enable_gates = False
                sampler.options.twirling.enable_measure = False

            protein_solver = ProteinSolver(
                ansatz=isa_circ,
                hamiltonian=qubit_op,
                sampler=sampler,
                parallelizer=PARALLELIZER,
            )

            opt_results = protein_solver.train(
                optimizer=OPTIMIZER,
                maxiter=MAX_ITER,
                num_batches=num_workers,
                init_params=init_params,
                num_saved_states=MAX_NUM_SAVED_STATES,
                verbose=True,
            )

    metadata["init_param_interval"] = "[4.71, 6.28]"
    metadata["run_id"] = run_id

    timestamp_str = TIMESTAMP.strftime("%Y-%m-%d-%H-%M-%S")
    out_dir = os.path.join(base_save_dir, f"{MAIN_SEQ}_{timestamp_str}_mps_run{run_id}")
    os.makedirs(out_dir, exist_ok=True)

    with open(os.path.join(out_dir, "metadata.json"), "w") as jf:
        json.dump(metadata, jf, indent=2, sort_keys=True)

    with open(os.path.join(out_dir, "opt_results.json"), "w") as jf:
        json.dump(opt_results, jf, indent=2)

    print(f"✅ Results saved to {out_dir}")

ray.shutdown()
