import sys
sys.path.append("../")

from qiskit.circuit.library import RealAmplitudes
from qiskit.quantum_info import SparsePauliOp
from qiskit_ibm_runtime import SamplerV2 as Sampler
from qiskit_ibm_runtime.fake_provider import FakeSherbrooke
from qiskit_ibm_runtime import QiskitRuntimeService, Session
from qiskit_aer import AerSimulator
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
import numpy as np
import json
import os
from datetime import datetime, timezone
from time import time

# ============================ 
# Define parameters 
MAIN_SEQ = "KLVFFA"
INIT_PARAMS_SOURCE = "KLVFFA_2025-05-24-16-44-24_aer_simulator_matrix_product_state_interval4"  
NUM_WORKERS = None
QUBIT_NUMBER = 24
ANSATZ_REPS = 2
RUNNER = "hardware"  # "aer-sv", "aer-mps", "hardware"
BACKEND_NAME = "ibm_cleveland"
SHOTS = 100_000
TIMESTAMP = datetime.now(timezone.utc)
# ============================

# Metadata tracking
metadata = {
    "main_seq": MAIN_SEQ,
    "ansatz_reps": ANSATZ_REPS,
    "runner": RUNNER,
    "shots": SHOTS,
    "timestamp": TIMESTAMP.isoformat(),
    "init_param_source": INIT_PARAMS_SOURCE,
}


def load_opt_params(folder_name: str) -> np.ndarray:
    path = f"../res/fcc_hw/{folder_name}/opt_results.json"
    if not os.path.exists(path):
        raise FileNotFoundError(f"Cannot find file: {path}")
    with open(path, "r") as f:
        data = json.load(f)
    params = data.get("opt_params", None)
    if params is None:
        raise ValueError(f"`opt_params` not found in {path}")
    return np.array(params)


if __name__ == "__main__":
    print(f"📥 Loading optimal parameters from: {INIT_PARAMS_SOURCE}")
    init_params = load_opt_params(INIT_PARAMS_SOURCE)

    print(f"⚙️ Building ansatz with {ANSATZ_REPS} repetitions...")
    ansatz = RealAmplitudes(QUBIT_NUMBER, reps=ANSATZ_REPS).decompose()
    ansatz.measure_all()
    
    print(f"🚀 Preparing backend: {RUNNER}")
    time_start = time()
    if RUNNER == "aer-sv":
        backend = AerSimulator(method="statevector")
    elif RUNNER == "aer-mps":
        backend = AerSimulator(method="matrix_product_state")
    elif RUNNER == "hardware":
        service = QiskitRuntimeService()
        backend = service.backend(BACKEND_NAME)
    else:
        raise ValueError(f"Unknown runner: {RUNNER}")
    metadata["backend"] = backend.name

    # Transpile circuit
    pass_manager = generate_preset_pass_manager(backend=backend, optimization_level=3)
    isa_circ = pass_manager.run(ansatz)
    time_end = time()
    metadata["circuit_prep_time (s)"] = round(time_end - time_start, 2)

    # Extract final physical qubit mapping (used qubit indices)
    physical_qubits_used = isa_circ.layout.final_index_layout()
    metadata["physical_qubits_used"] = physical_qubits_used


    # Run sampling using Qiskit Runtime
    with Session(backend=backend) as session:
        print("📡 Starting Qiskit Runtime session...")
        sampler = Sampler(mode=session)
        sampler.options.default_shots = SHOTS

        if RUNNER == "hardware":
            sampler.options.twirling.enable_gates = True
            sampler.options.twirling.enable_measure = True

        print("🎯 Submitting sampling job...")
        job = sampler.run([(isa_circ, init_params)])
        result = job.result()
        metadata["session_id"] = session.session_id
        print("✅ Sampling complete!")

        # Save results and metadata
        timestamp_str = TIMESTAMP.strftime("%Y-%m-%d-%H-%M-%S")
        out_dir = f"../res/fcc_hw/{MAIN_SEQ}_{timestamp_str}_{backend.name}_sampled"
        os.makedirs(out_dir, exist_ok=True)

        with open(f"{out_dir}/sampling_results.json", "w") as f:
            json.dump(result[0].data.meas.get_counts(), f, indent=2)

        with open(f"{out_dir}/metadata.json", "w") as f:
            json.dump(metadata, f, indent=2, sort_keys=True)

        print(f"📁 Results saved to {out_dir}")
