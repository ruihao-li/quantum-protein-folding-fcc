from __future__ import annotations
from vqec import VQECResult
from qiskit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp
from qiskit.primitives import BaseEstimator
from qiskit_algorithms.gradients import BaseEstimatorGradient
import numpy as np
from tqdm import tqdm

class VQECOpt_goldstein:
    def __init__(
        self,
        qubit_op: SparsePauliOp,
        constr_ops: list[SparsePauliOp],
        ansatz: QuantumCircuit,
        estimator: BaseEstimator,
        gradient: BaseEstimatorGradient,
        initial_params: np.ndarray | None = None,
        initial_dual_vars: np.ndarray | None = None,
        primal_perturb_step: float = 0.05,
        dual_perturb_step: float = 0.05,
        gamma: float = 1.0,
        max_iter: int = 200,
        alpha_min: float = 1e-6,   # Smallest allowed primal step size
        alpha_max: float = 1.0,    # Maximum allowed primal step size
        omega: float = 0.2,        # Goldstein test tolerance
        theta: float = 0.5         # Reduction factor for step size
    ):
        self._qubit_op = qubit_op
        self._constr_ops = constr_ops
        self._ansatz = ansatz
        self._estimator = estimator
        self._gradient = gradient
        self._initial_params = (
            initial_params if initial_params is not None else np.random.uniform(0, 2 * np.pi, ansatz.num_parameters)
        )
        self._initial_dual_vars = (
            initial_dual_vars if initial_dual_vars is not None else np.array([0.1] * len(constr_ops))
        )
        self._primal_perturb_step = primal_perturb_step
        self._dual_perturb_step = dual_perturb_step
        self._gamma = gamma
        self._max_iter = max_iter
        self._alpha_min = alpha_min
        self._alpha_max = alpha_max
        self._omega = omega
        self._theta = theta

    def get_expectation(self, circuit: QuantumCircuit, observable: SparsePauliOp, params: np.ndarray) -> float:
        return self._estimator.run(circuit, observable, params).result().values[0]

    def perturb_primal(self, x, y, obj_grad, constr_grads):
        # Perturb primal variables according to equation (1.4)
        primal_var_perturbations = (y * constr_grads.T).T.sum(axis=0) + obj_grad
        perturbed_primal_vars = x - self._primal_perturb_step * primal_var_perturbations
        return np.mod(perturbed_primal_vars, 4 * np.pi)

    def perturb_dual(self, x, y, constraints):
        # Perturb dual variables according to equation (1.5)
        perturbed_dual_vars = np.maximum(y + self._dual_perturb_step * constraints, 0)
        return perturbed_dual_vars

    def evaluate_gap_function(self, primal_vars, dual_vars, perturbed_primal_vars):
        # Compute Lagrangian values and gap function (Eq. 1.6 and 1.7)
        lagrangian_perturbed_dual = self.get_expectation(self._ansatz, self._qubit_op, primal_vars) + np.dot(dual_vars, self.get_constraints(primal_vars))
        lagrangian_perturbed_primal = self.get_expectation(self._ansatz, self._qubit_op, perturbed_primal_vars) + np.dot(dual_vars, self.get_constraints(perturbed_primal_vars))
        gap_primal = lagrangian_perturbed_dual - lagrangian_perturbed_primal
        return gap_primal

    def get_constraints(self, primal_vars):
        return np.array([self.get_expectation(self._ansatz, constr_op, primal_vars) for constr_op in self._constr_ops])

    def apply_goldstein_test(self, primal_vars, dual_vars, perturbed_primal_vars):
        # Calculate the derivative and perform the Goldstein test
        gap_primal = self.evaluate_gap_function(primal_vars, dual_vars, perturbed_primal_vars)
        if gap_primal >= self._omega * self._primal_perturb_step * np.linalg.norm(perturbed_primal_vars - primal_vars):
            return True  # Test passes
        return False

    def optimize_primal_dual(self):
        primal_vars = np.array([self._initial_params])
        dual_vars = self._initial_dual_vars

        energies = [self.get_expectation(self._ansatz, self._qubit_op, primal_vars)]
        constraints = [self.get_constraints(primal_vars)]
        lagrangian_gap_trajectory = []

        for i in tqdm(range(self._max_iter)):
            # Compute gradients for objective and constraints
            obj_grad = np.array(self._gradient.run(circuits=self._ansatz, observables=self._qubit_op, parameter_values=primal_vars).result().gradients)
            constr_grads = np.array([self._gradient.run(circuits=self._ansatz, observables=constr_op, parameter_values=primal_vars).result().gradients for constr_op in self._constr_ops])

            perturbed_primal_vars = self.perturb_primal(primal_vars, dual_vars, obj_grad, constr_grads)
            perturbed_dual_vars = self.perturb_dual(primal_vars, dual_vars, constraints[-1])

            if not self.apply_goldstein_test(primal_vars, dual_vars, perturbed_primal_vars):
                # Reduce step size and retry perturbation
                print(f"Goldstein test failed at iteration {i}. Reducing step size.")
                self._primal_perturb_step = max(self._alpha_min, (1 - self._theta) * self._primal_perturb_step)
                continue  # Go back and retry with new step size

            # If test passes, update primal and dual variables
            lagrangian_gap = self.evaluate_gap_function(primal_vars, dual_vars, perturbed_primal_vars)
            lagrangian_gap_trajectory.append(lagrangian_gap)

            grad_primal = -(obj_grad + (perturbed_dual_vars * constr_grads.T).T.sum(axis=0))
            grad_dual = self.get_constraints(perturbed_primal_vars)

            step_size = self._gamma * lagrangian_gap / (np.linalg.norm(grad_primal) ** 2 + np.linalg.norm(grad_dual) ** 2)
            primal_vars = primal_vars + step_size * grad_primal
            primal_vars = np.mod(primal_vars, 4 * np.pi)

            dual_vars = np.maximum(dual_vars + step_size * grad_dual, 0)

            energies.append(self.get_expectation(self._ansatz, self._qubit_op, primal_vars))
            constraints.append(self.get_constraints(primal_vars))

        return self._build_vqec_result(energies, constraints, primal_vars, dual_vars, lagrangian_gap_trajectory)

    def _build_vqec_result(self, energies, constraints, primal_vars, dual_vars, lagrangian_gaps):
        result = VQECResult()
        result.ansatz = self._ansatz.copy()
        result.energies = np.array(energies)
        result.constraints = np.array(constraints)
        result.optimal_primal_vars = primal_vars[0]
        result.optimal_dual_vars = dual_vars
        result.vqec_iterations = len(energies) - 1
        result.lagrangian_gap_trajectory = np.array(lagrangian_gaps)
        return result
