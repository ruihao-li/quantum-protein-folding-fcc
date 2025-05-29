import sys

sys.path.append("../")

from qiskit.circuit.library import RealAmplitudes, TwoLocal
from qiskit.quantum_info import SparsePauliOp
from qiskit_ibm_runtime import SamplerV2 as Sampler
from qiskit_ibm_runtime.fake_provider import FakeSherbrooke
from qiskit_ibm_runtime import QiskitRuntimeService, Session
from qiskit_aer import AerSimulator
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
import collections
import numpy as np
import json
import os
from datetime import datetime, timezone
from time import time
import psutil

# ============================
# Define parameters
MAIN_SEQ = "KLVFFA"
OPT_PARAMS_SOURCE = "klvffa_vqec_p0_5_u2/klvffa_vqec_res_12.json"
QUBIT_NUMBER = 24
ANSATZ_REPS = 2
RUNNER = "hardware"  # "aer-sv", "aer-mps", "hardware"
BACKEND_NAME = "ibm_cleveland"  # change to other backends if needed, comment out if not using hardware
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
    "vqec_opt_params_source": OPT_PARAMS_SOURCE,
}


def load_opt_params(source_path: str) -> np.ndarray:
    """Load optimal parameters from the simulation run."""
    path = f"../notebooks/vqec_results/{source_path}"
    if not os.path.exists(path):
        raise FileNotFoundError(f"Cannot find file: {path}")
    with open(path, "r") as f:
        data = json.load(f)
    params = data.get("optimal_primal_vars", None)
    if params is None:
        raise ValueError(f"`optimal_primal_vars` not found in {path}")
    return np.array(params)


if __name__ == "__main__":
    # Load parameters
    print(f"📥 Loading optimal parameters from: {OPT_PARAMS_SOURCE}")
    opt_params = load_opt_params(OPT_PARAMS_SOURCE)

    # Create ansatz circuit
    print(f"⚙️ Building ansatz with {ANSATZ_REPS} repetitions...")
    ansatz = TwoLocal(
        QUBIT_NUMBER,
        ["ry"],
        "cx",
        reps=ANSATZ_REPS,
        entanglement="linear",
        skip_final_rotation_layer=True,
    )
    ansatz.measure_all()

    # Backend selection
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

    # Run sampling using Qiskit Runtime
    with Session(backend=backend) as session:
        print("📡 Starting Qiskit Runtime session...")
        sampler = Sampler(mode=session)
        sampler.options.default_shots = SHOTS

        if RUNNER == "hardware":
            sampler.options.twirling.enable_gates = True
            sampler.options.twirling.enable_measure = True

        print("🎯 Submitting sampling job...")
        job = sampler.run([(isa_circ, opt_params)])
        result = job.result()
        # pub_result = result[0].data
        # counts = pub_result.meas.get_counts()
        metadata["session_id"] = session.session_id
        print("✅ Sampling complete!")

        # Save sampler results and metadata
        timestamp_str = TIMESTAMP.strftime("%Y-%m-%d-%H-%M-%S")
        out_dir = f"../res/vqec_hw/{MAIN_SEQ}_{timestamp_str}_{backend.name}_sampled"
        os.makedirs(out_dir, exist_ok=True)

        with open(f"{out_dir}/sampling_results.json", "w") as f:
            json.dump(result[0].data.meas.get_counts(), f, indent=2)

        with open(f"{out_dir}/metadata.json", "w") as f:
            json.dump(metadata, f, indent=2, sort_keys=True)

        print(f"📁 Results saved to {out_dir}")
