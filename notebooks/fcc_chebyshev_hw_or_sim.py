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
MAIN_SEQ = "GNLVS"
NUM_WORKERS = None
R2_THRESHOLD = 0.999
CHUNK = 20
RUNNER = "aer"  # "aer" or "hardware"
PARALLELIZER = "ray"  # "ray" or "python-mp"
SHOTS = 50_000
OPTIMIZER = "COBYLA"
MAX_ITER = 10
MAX_NUM_SAVED_STATES = 1000
TIMESTAMP = datetime.now(timezone.utc)
# ============================

metadata["main_seq"] = MAIN_SEQ
metadata["r2_threshold"] = R2_THRESHOLD
metadata["chunk size"] = CHUNK
metadata["runner"] = RUNNER
metadata["parallelizer"] = PARALLELIZER
metadata["shots"] = SHOTS
metadata["optimizer"] = OPTIMIZER
metadata["max_iter"] = MAX_ITER
metadata["max_num_saved_states"] = MAX_NUM_SAVED_STATES
metadata["timestamp"] = TIMESTAMP.isoformat()


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


if __name__ == "__main__":

    ####################
    num_workers = (
        NUM_WORKERS if NUM_WORKERS is not None else psutil.cpu_count(logical=False)
    )
    print(f"Number of workers: {num_workers}")
    metadata["num_workers"] = num_workers

    if PARALLELIZER == "ray":
        ray.init(
            num_cpus=num_workers,
            ignore_reinit_error=True,
            log_to_driver=False,
            runtime_env={
                "py_modules": [fcc],
            },
        )

    ham_dir = "fcc_hamilts"
    try:
        os.makedirs(ham_dir)
    except FileExistsError:
        pass
    # Check if the Hamiltonian operator is already generated
    if os.path.exists(f"{ham_dir}/{MAIN_SEQ}_{len(MAIN_SEQ)}AAs.json"):
        print("Hamiltonian operator found.")
        tic = time()
        with open(f"{ham_dir}/{MAIN_SEQ}_{len(MAIN_SEQ)}AAs.json", "r") as f:
            data = json.load(f)
            num_qubits = data["num_qubits"]
            sparse_op_list = [
                (pauli, qubits, coeff)
                for pauli, qubits, coeff in data["sparse_op_list"]
            ]

        qubit_op = SparsePauliOp.from_sparse_list(sparse_op_list, num_qubits=num_qubits)
        toc = time()
        ham_load_time = toc - tic  # in seconds
        print(f"Hamiltonian loaded in {ham_load_time / 60:.2f} min")
        metadata["hamiltonian_load_time (min)"] = np.round(ham_load_time / 60, 2)
    else:
        print("Hamiltonian operator not found. Generating...")
        tic = time()
        pf_problem = build_pf(MAIN_SEQ)
        qubit_op = pf_problem.qubit_op(r2_threshold=R2_THRESHOLD, chunk=CHUNK)
        toc = time()
        ham_gen_time = toc - tic  # in seconds
        print(f"Hamiltonian generated in {ham_gen_time / 60:.2f} min")
        metadata["hamiltonian_gen_time (min)"] = np.round(ham_gen_time / 60, 2)
        # Save the Hamiltonian operator to a file
        with open(f"fcc_hamilts/{MAIN_SEQ}_{len(MAIN_SEQ)}AAs.json", "w") as f:
            # Convert complex coefficients to real numbers
            qubit_op_real = qubit_op.copy()
            qubit_op_real.coeffs.dtype = np.float64
            sparse_op_list = qubit_op_real.to_sparse_list()
            num_qubits = qubit_op_real.num_qubits
            json.dump(
                {"num_qubits": num_qubits, "sparse_op_list": sparse_op_list},
                f,
            )
        print("Hamiltonian operator saved.")
    print(f"Number of qubits: {qubit_op.num_qubits}")
    print(f"Number of Hamiltonian terms: {len(qubit_op)}")

    metadata["num_qubits"] = qubit_op.num_qubits
    metadata["num_ham_terms"] = len(qubit_op)

    ansatz = RealAmplitudes(qubit_op.num_qubits, reps=1).decompose()
    # Add measurements
    ansatz.measure_all()
    # Get backend
    time_start = time()

    if RUNNER == "aer":
        backend = AerSimulator(
            method="matrix_product_state"
        )  # MPS by default (manually change it to "statevector" for smaller problems)
    elif RUNNER == "hardware":
        service = QiskitRuntimeService()
        backend = service.backend("ibm_cleveland")
        # backend = FakeSherbrooke()
    print(f"Backend: {backend.name}")
    metadata["backend"] = backend.name
    # Transpilation for hardware runs
    pass_manager = generate_preset_pass_manager(backend=backend, optimization_level=3)
    isa_circ = pass_manager.run(ansatz)
    time_end = time()
    print(f"Time to prepare the circuit: {time_end - time_start:.2f}s")
    metadata["circ_preparation_time (s)"] = np.round(time_end - time_start, 2)

    # Run the VQE

    with Session(backend=backend) as session:
        metadata["session_id"] = session.session_id
        sampler = Sampler(mode=session)
        sampler.options.default_shots = SHOTS
        if RUNNER == "hardware":
            # recommended: for Sampler runs enable gates twirling and disable measure twirling
            sampler.options.twirling.enable_gates = False
            sampler.options.twirling.enable_measure = False
        print("Starting VQE with parallelizer:", PARALLELIZER)
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

# ray.shutdown()
