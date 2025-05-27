import sys

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
from qiskit.quantum_info import SparsePauliOp
from qiskit_ibm_runtime import SamplerV2 as Sampler
from qiskit_ibm_runtime.fake_provider import FakeSherbrooke  # for local testing mode
from qiskit_ibm_runtime import QiskitRuntimeService, Session
from qiskit_aer import AerSimulator
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from math import pi
import ray
import psutil
import json
from time import time
import numpy as np
from datetime import datetime, timezone
import os

metadata = {}
# ============================
# Define parameters
MAIN_SEQ = "KLVFFA"
NUM_WORKERS = None
R2_THRESHOLD = 0.999
CHUNK = None
ANSATZ_REPS = 2
INIT_PARAMS_FILE = None
RUNNER = "hardware"
PARALLELIZER = "ray"
SHOTS = 100_000
OPTIMIZER = "COBYLA"
MAX_ITER = 500
MAX_NUM_SAVED_STATES = 1000
TIMESTAMP = datetime.now(timezone.utc)
# ============================

metadata["main_seq"] = MAIN_SEQ
metadata["r2_threshold"] = R2_THRESHOLD
metadata["chunk_size"] = CHUNK
metadata["ansatz_reps"] = ANSATZ_REPS
metadata["init_params_file"] = INIT_PARAMS_FILE
metadata["runner"] = RUNNER
metadata["parallelizer"] = PARALLELIZER
metadata["shots"] = SHOTS
metadata["optimizer"] = OPTIMIZER
metadata["max_iter"] = MAX_ITER
metadata["max_num_saved_states"] = MAX_NUM_SAVED_STATES
metadata["timestamp"] = TIMESTAMP.isoformat()


def build_pf(main_seq: str, energy_matrix_file: str = "mj_matrix"):
    mj_interaction = MiyazawaJerniganInteraction(energy_matrix_file)
    penalty_terms = PenaltyParameters(50, 50, 50)
    peptide = Peptide(main_seq)
    return ProteinFoldingProblem(peptide, mj_interaction, penalty_terms)

def load_init_params(init_params_file: str | None) -> np.ndarray | None:
    """
    Load initial parameters from the `opt_results` file (typically from an MPS
    run).

    Args:
        init_params_file (str): Path to the `opt_results` file.

    Returns:
        list: Initial parameters to be used in the VQE optimization.
    """ 
    if init_params_file is None:
        return None
    init_params_file_path = f"../res/fcc_hw/{init_params_file}/opt_results.json"
    if not os.path.exists(init_params_file_path):
        raise FileNotFoundError(f"File not found: {init_params_file_path}")
    with open(init_params_file_path, "r") as f:
        data = json.load(f)
    init_params = data.get("init_params", None)
    init_params = np.array(init_params) if init_params is not None else None
    return init_params

if __name__ == "__main__":
    num_workers = NUM_WORKERS if NUM_WORKERS is not None else psutil.cpu_count(logical=False)
    print(f"Number of workers: {num_workers}")
    metadata["num_workers"] = num_workers

    if PARALLELIZER == "ray":
        ray.init(
            num_cpus=num_workers,
            ignore_reinit_error=True,
            log_to_driver=False,
            runtime_env={"py_modules": [fcc]},
        )

    ham_dir = "fcc_hamilts"
    file_name = f"{MAIN_SEQ}_{len(MAIN_SEQ)}AAs_R2_{int(R2_THRESHOLD * 1000)}.json"
    os.makedirs(ham_dir, exist_ok=True)

    if os.path.exists(f"{ham_dir}/{file_name}"):
        print("Hamiltonian operator found.")
        tic = time()
        with open(f"{ham_dir}/{file_name}", "r") as f:
            data = json.load(f)
            num_qubits = data["num_qubits"]
            sparse_op_list = [
                (label, complex(coeff["real"], coeff["imag"]))
                for label, coeff in data["sparse_list"]
            ]
        qubit_op = SparsePauliOp.from_list(sparse_op_list, num_qubits=num_qubits)
        toc = time()
        metadata["hamiltonian_load_time (min)"] = np.round((toc - tic) / 60, 2)
    else:
        print("Hamiltonian operator not found. Generating...")
        tic = time()
        pf_problem = build_pf(MAIN_SEQ)
        qubit_op = pf_problem.qubit_op(r2_threshold=R2_THRESHOLD, chunk=CHUNK)
        toc = time()
        metadata["hamiltonian_gen_time (min)"] = np.round((toc - tic) / 60, 2)
        with open(f"{ham_dir}/{file_name}", "w") as f:
            qubit_op_real = qubit_op.copy()
            qubit_op_real.coeffs.dtype = np.float64
            sparse_op_list = qubit_op_real.to_sparse_list()
            json.dump(
                {"num_qubits": qubit_op_real.num_qubits, "sparse_list": sparse_op_list},
                f,
            )

    print(f"Number of qubits: {qubit_op.num_qubits}")
    print(f"Number of Hamiltonian terms: {len(qubit_op)}")
    metadata["num_qubits"] = qubit_op.num_qubits
    metadata["num_ham_terms"] = len(qubit_op)

    ansatz = RealAmplitudes(qubit_op.num_qubits, reps=ANSATZ_REPS).decompose()
    ansatz.measure_all()

    time_start = time()
    if RUNNER == "aer-sv":
        backend = AerSimulator(method="statevector")
    elif RUNNER == "aer-mps":
        backend = AerSimulator(method="matrix_product_state")
    elif RUNNER == "hardware":
        service = QiskitRuntimeService()
        backend = service.backend("ibm_cleveland")
    print(f"Backend: {backend.name}")
    metadata["backend"] = backend.name

    pass_manager = generate_preset_pass_manager(backend=backend, optimization_level=3)
    isa_circ = pass_manager.run(ansatz)
    time_end = time()
    metadata["circ_preparation_time (s)"] = np.round(time_end - time_start, 2)

    # === Initial parameters sampled from interval 2: [π/2, π] ===
    num_params = isa_circ.num_parameters
    init_params = load_init_params("KLVFFA_2025-05-15-19-51-03_aer_simulator_matrix_product_state_interval2")

    with Session(backend=backend) as session:
        metadata["session_id"] = session.session_id
        sampler = Sampler(mode=session)
        sampler.options.default_shots = SHOTS
        if RUNNER == "hardware":
            sampler.options.twirling.enable_gates = True
            sampler.options.twirling.enable_measure = True

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

        metadata["init_param_interval"] = "[π/2, π]"
        metadata["interval_index"] = 2

        timestamp_str = TIMESTAMP.strftime("%Y-%m-%d-%H-%M-%S")
        parent_dir = f"../res/fcc_hw/{MAIN_SEQ}_{timestamp_str}_{backend.name}_interval2_hw"
        os.makedirs(parent_dir, exist_ok=True)

        with open(f"{parent_dir}/metadata.json", "w") as jf:
            json.dump(metadata, jf, indent=2, sort_keys=True)

        with open(f"{parent_dir}/opt_results.json", "w") as jf:
            json.dump(opt_results, jf, indent=2)

        print(f"✅ Results saved to {parent_dir}")

    ray.shutdown()
