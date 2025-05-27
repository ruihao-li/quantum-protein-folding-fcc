import sys

sys.path.append("..")
from fcc import (
    MiyazawaJerniganInteraction,
    Peptide,
    ProteinFoldingProblem,
    PenaltyParameters,
)
from vqe import PerturbedPrimalDualOpt, VQECResult
from qiskit.circuit.library import RealAmplitudes, EfficientSU2, TwoLocal

# from qiskit_aer.primitives import Estimator, Sampler
from qiskit_aer.primitives import Estimator
from qiskit_algorithms.gradients import (
    LinCombEstimatorGradient,
    ParamShiftEstimatorGradient,
    FiniteDiffEstimatorGradient,
)
import numpy as np
import ray
import psutil
import os
from time import time
import vqe


# ============================
# Define parameters
MAIN_SEQ = "KLVFFA"
NUM_WORKERS = 30  # 48
ANSATZ_REPS = 2
MAX_ITER = 1500  # 500
NUM_OPT_REPS = 20  # 20
PRIMAL_PERTURB_STEP = 0.2
DUAL_PERTURB_STEP = 0.2
PRIMAL_DUAL_UPDATE_STEP = 0.1

# Grid search perturbation and update steps
# PERTURB_STEPS_LIST = [0.01, 0.05, 0.1, 0.2, 0.5]
# UPDATE_STEPS_LIST = [0.1, 0.5, 1, 2, 5]
# ============================


def build_pf(main_seq: str, energy_matrix_file: str = "mj_matrix"):
    """Builds the protein folding problem for the given sequence."""

    mj_interaction = MiyazawaJerniganInteraction(energy_matrix_file)
    print(mj_interaction.calculate_energy_matrix(main_seq))

    penalty_back = 50
    penalty_redun = 50
    penalty_olap = None

    penalty_terms = PenaltyParameters(
        penalty_back=penalty_back,
        penalty_redun=penalty_redun,
        penalty_olap=penalty_olap,
    )

    peptide = Peptide(main_seq)

    protein_folding_problem = ProteinFoldingProblem(
        peptide, mj_interaction, penalty_terms
    )

    return protein_folding_problem


@ray.remote
def _optimize_ray(
    ppd_opt: PerturbedPrimalDualOpt,
    init_params: np.ndarray | None = None,
    init_dual_vars: np.ndarray | None = None,
    primal_perturb_step: float = 0.05,
    dual_perturb_step: float = 0.05,
    primal_dual_update_step: float = 1,
    auto_update_step: bool = False,
    max_iter: int = 300,
) -> VQECResult:
    """Optimize the given problem using the Perturbed Primal Dual method."""
    return ppd_opt.optimize_primal_dual(
        initial_params=init_params,
        initial_dual_vars=init_dual_vars,
        primal_perturb_step=primal_perturb_step,
        dual_perturb_step=dual_perturb_step,
        gamma=primal_dual_update_step,
        auto_update_step=auto_update_step,
        max_iter=max_iter,
    )


# Initialize Ray
# num_cpus = psutil.cpu_count(logical=False)
# print(f"Number of CPUs: {num_cpus}")
ray.init(
    ignore_reinit_error=True,
    num_cpus=NUM_WORKERS,
    runtime_env={
        "py_modules": [vqe],
    },
)

protein_folding_problem = build_pf(MAIN_SEQ)
qubit_op = protein_folding_problem.qubit_op()
print(f"Number of terms: {len(qubit_op)}")
print(f"Number of qubits: {qubit_op.num_qubits}")

bead_pairs, olap_constr_ops = protein_folding_problem.olap_constr_ops()
print(f"Number of overlap constraints: {len(olap_constr_ops)}")
print(f"Number of qubits: {olap_constr_ops[0].num_qubits}")

ansatz = TwoLocal(
    qubit_op.num_qubits,
    ["ry"],
    "cx",
    reps=ANSATZ_REPS,
    entanglement="linear",
    skip_final_rotation_layer=True,
)
estimator = Estimator(
    backend_options={"method": "matrix_product_state"},
    run_options={"shots": None},
    approximation=True,
)
gradient = ParamShiftEstimatorGradient(estimator=estimator)

start_time = time()
init_dual_vars = np.array([0.0] * len(olap_constr_ops))  # Fix initial dual vars
vqec = PerturbedPrimalDualOpt(
    qubit_op=qubit_op,
    constr_ops=olap_constr_ops,
    ansatz=ansatz,
    estimator=estimator,
    gradient=gradient,
)

res = ray.get(
    [
        _optimize_ray.remote(
            ppd_opt=vqec,
            init_params=None,
            init_dual_vars=init_dual_vars,
            max_iter=MAX_ITER,
            auto_update_step=False,
            primal_perturb_step=PRIMAL_PERTURB_STEP,
            dual_perturb_step=DUAL_PERTURB_STEP,
            primal_dual_update_step=PRIMAL_DUAL_UPDATE_STEP,
        )
        for _ in range(NUM_OPT_REPS)
    ]
)

perturb_step = str(PRIMAL_PERTURB_STEP).replace(".", "_")
update_step = str(PRIMAL_DUAL_UPDATE_STEP).replace(".", "_")
directory = f"klvffa_vqec_p{perturb_step}_u{update_step}"
os.makedirs("vqec_results/" + directory, exist_ok=True)
# Save results to JSON files
for i, r in enumerate(res):
    if i < 10:
        file_name = f"{directory}/klvffa_vqec_res_0{i}.json"
    else:
        file_name = f"{directory}/klvffa_vqec_res_{i}.json"
    r.write_to_json(file_name)
    print(f"Result {i} saved to {file_name}")

end_time = time()
print(f"Total time taken: {(end_time - start_time) / 60:.2f} minutes")

# Shut down Ray
ray.shutdown()
