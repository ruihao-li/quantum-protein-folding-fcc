import sys
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
from qiskit_ibm_runtime import SamplerV2 as Sampler
from qiskit_ibm_runtime import QiskitRuntimeService, Session
from qiskit_aer import AerSimulator
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from math import pi
import ray
import psutil
import json
import numpy as np
from datetime import datetime, timezone
import os

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

# === Helper Function ===
def build_pf(main_seq: str, energy_matrix_file: str = "mj_matrix"):
    mj_interaction = MiyazawaJerniganInteraction(energy_matrix_file)
    penalty_terms = PenaltyParameters(50, 50, 50)
    peptide = Peptide(main_seq)
    return ProteinFoldingProblem(peptide, mj_interaction, penalty_terms)

# === Main Execution ===
if __name__ == "__main__":
    num_workers = NUM_WORKERS or psutil.cpu_count(logical=False)
    metadata["num_workers"] = num_workers
    if PARALLELIZER == "ray":
        ray.init(
            num_cpus=num_workers,
            ignore_reinit_error=True,
            log_to_driver=False,
            runtime_env={"py_modules": [fcc]},
        )

    # Load or generate Hamiltonian
    ham_dir = "fcc_hamilts"
    file_name = f"{MAIN_SEQ}_{len(MAIN_SEQ)}AAs_R2_{int(R2_THRESHOLD * 1000)}.json"
    os.makedirs(ham_dir, exist_ok=True)

    if os.path.exists(f"{ham_dir}/{file_name}"):
        with open(f"{ham_dir}/{file_name}", "r") as f:
            data = json.load(f)
            num_qubits = data["num_qubits"]
            sparse_op_list = [(label, coeff) for label, coeff in data["sparse_list"]]
        qubit_op = SparsePauliOp.from_list(sparse_op_list, num_qubits=num_qubits)
    else:
        pf_problem = build_pf(MAIN_SEQ)
        qubit_op = pf_problem.qubit_op(r2_threshold=R2_THRESHOLD, chunk=CHUNK)
        with open(f"{ham_dir}/{file_name}", "w") as f:
            qubit_op_real = qubit_op.copy()
            qubit_op_real.coeffs.dtype = np.float64
            sparse_op_list = qubit_op_real.to_sparse_list()
            json.dump({"num_qubits": qubit_op_real.num_qubits, "sparse_list": sparse_op_list}, f)

    metadata["num_qubits"] = qubit_op.num_qubits
    metadata["num_ham_terms"] = len(qubit_op)

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

    intervals = [
        (pi, 3 * pi / 2),
        (3 * pi / 2, 2 * pi),
    ]

    for idx, (low, high) in enumerate(intervals, 3):
        print(f"\n--- Running VQE for Interval {idx}: [{low:.2f}, {high:.2f}] ---")
        init_params = np.random.uniform(low, high, size=isa_circ.num_parameters)

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

            metadata["init_param_interval"] = f"[{low:.2f}, {high:.2f}]"
            metadata["interval_index"] = idx

            timestamp_str = TIMESTAMP.strftime("%Y-%m-%d-%H-%M-%S")
            out_dir = f"../res/fcc_hw/{MAIN_SEQ}_{timestamp_str}_{backend.name}_interval{idx}"
            os.makedirs(out_dir, exist_ok=True)

            with open(f"{out_dir}/metadata.json", "w") as jf:
                json.dump(metadata, jf, indent=2, sort_keys=True)

            with open(f"{out_dir}/opt_results.json", "w") as jf:
                json.dump(opt_results, jf, indent=2)

            print(f"✅ Results saved to {out_dir}")

    ray.shutdown()
