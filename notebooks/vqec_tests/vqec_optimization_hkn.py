"""Variational Quantum Eigensolver with Constraints (VQEC) module. Based on
arXiv:2311.08502. This version contains an automatic selection of the update
step size based on the method proposed in: M. Kallio and A. Ruszczynski,
Perturbation methods for saddle point computation, Tech. Rep. (International
Institute for Applied Systems Analysis, Laxenburg, Austria: WP-94-038, 1994)."""

from __future__ import annotations
from vqec import VQECResult
from qiskit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp
from qiskit.primitives import BaseEstimator
from qiskit_algorithms.gradients import BaseEstimatorGradient
import numpy as np
from tqdm import tqdm


class VQECOpt_hkn:

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
    ):
        """
        Args:
            qubit_op: The qubit operator for which the expectation value is minimized.
            constr_ops: The constraint operators whose the expectation values are constrained.
            ansatz: The quantum circuit that prepares the quantum state.
            estimator: The estimator that computes the expectation values.
            gradient: The gradient tool that computes the gradients of the expectation values.
            initial_params: The initial parameters of the quantum circuit.
            initial_dual_vars: The initial dual variables of the optimization.
            primal_perturb_step: The perturbation step size for the primal variables.
            dual_perturb_step: The perturbation step size for the dual variables.
            gamma: The step size parameter used for the primal-dual update.
            max_iter: The maximum number of iterations of the optimization.
        """
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

    def get_expectation(
        self, circuit: QuantumCircuit, observable: SparsePauliOp, params: np.ndarray
    ) -> float:
        """
        Computes the expectation value of an observable with respect to a quantum state.

        Args:
            circuit: The quantum circuit that prepares the quantum state.
            observable: The observable for which the expectation value is computed.
            params: The parameters of the quantum circuit.

        Returns:
            The expectation value of the observable.
        """
        return self._estimator.run(circuit, observable, params).result().values[0]

    def optimize_primal_dual(self) -> VQECResult:
        """
        Optimizes the primal and dual problems of the VQEC.

        Returns:
            The result of the optimization.
        """
        unperturbed_primal_vars = np.array([self._initial_params])
        unperturbed_dual_vars = self._initial_dual_vars

        energies = [
            self.get_expectation(self._ansatz, self._qubit_op, unperturbed_primal_vars)
        ]
        constraints = [
            [
                self.get_expectation(self._ansatz, constr_op, unperturbed_primal_vars)
                for constr_op in self._constr_ops
            ]
        ]
        lagrangian_gap_trajectory = []
        gradient_norms = []  # List to store the gradient norms

        for i in tqdm(range(self._max_iter)):
            # Compute gradients for the objective and constraint operators
            obj_grad = np.array(
                (
                    self._gradient.run(
                        circuits=self._ansatz,
                        observables=self._qubit_op,
                        parameter_values=unperturbed_primal_vars,
                    ).result().gradients
                )
            )
            constr_grads = np.array(
                [
                    self._gradient.run(
                        circuits=self._ansatz,
                        observables=constr_op,
                        parameter_values=unperturbed_primal_vars,
                    ).result().gradients
                    for constr_op in self._constr_ops
                ]
            )

            # Compute the perturbed primal variables
            primal_var_perturbations = (unperturbed_dual_vars * constr_grads.T).T.sum(axis=0) + obj_grad
            perturbed_primal_vars = (
                unperturbed_primal_vars - self._primal_perturb_step * primal_var_perturbations
            )
            perturbed_primal_vars = np.mod(perturbed_primal_vars, 4 * np.pi)

            # Compute the perturbed dual variables
            dual_var_perturbations = np.array(constraints[-1])
            perturbed_dual_vars = np.maximum(
                unperturbed_dual_vars + self._dual_perturb_step * dual_var_perturbations,
                0,
            )

            # Compute the Lagrangian at perturbed dual variables
            lagrangian_perturbed_dual = (
                energies[-1] + perturbed_dual_vars @ dual_var_perturbations
            )
            obj_exp_perturbed_primal = self.get_expectation(
                self._ansatz, self._qubit_op, perturbed_primal_vars
            )
            constr_exp_perturbed_primal = np.array(
                [
                    self.get_expectation(self._ansatz, constr_op, perturbed_primal_vars)
                    for constr_op in self._constr_ops
                ]
            )
            lagrangian_perturbed_primal = (
                obj_exp_perturbed_primal + unperturbed_dual_vars @ constr_exp_perturbed_primal
            )
            lagrangian_gap = lagrangian_perturbed_dual - lagrangian_perturbed_primal
            lagrangian_gap_trajectory.append(lagrangian_gap)
            if np.abs(lagrangian_gap) < 1e-10:
                print(f"Gap closed after {i+1} iterations; VQEC converged.")
                break

            # Compute gradients
            grad_primal = -(obj_grad + (perturbed_dual_vars * constr_grads.T).T.sum(axis=0))
            grad_dual = constr_exp_perturbed_primal

            # Calculate and store the gradient norm
            gradient_norm = np.sqrt(np.linalg.norm(grad_primal) ** 2 + np.linalg.norm(grad_dual) ** 2)
            gradient_norms.append(gradient_norm)

            # Update step size
            step_size = self._gamma * lagrangian_gap / (np.linalg.norm(grad_primal) ** 2 + np.linalg.norm(grad_dual) ** 2)

            # Update primal and dual variables
            updated_primal_vars = unperturbed_primal_vars + step_size * grad_primal
            updated_primal_vars = np.mod(updated_primal_vars, 4 * np.pi)

            updated_dual_vars = np.maximum(
                unperturbed_dual_vars + step_size * grad_dual, 0
            )

            energies.append(
                self.get_expectation(self._ansatz, self._qubit_op, updated_primal_vars)
            )
            constraints.append(
                [
                    self.get_expectation(self._ansatz, constr_op, updated_primal_vars)
                    for constr_op in self._constr_ops
                ]
            )

            # Check for convergence
            tol = 1e-5
            if np.abs(energies[-1] - energies[-2]) / np.abs(energies[-2]) < tol:
                print(f"Converged after {i+1} iterations")
                break

            # Update the primal and dual variables
            unperturbed_primal_vars = updated_primal_vars
            unperturbed_dual_vars = updated_dual_vars

        return self._build_vqec_result(
            energies=energies,
            constraints=constraints,
            primal_vars=unperturbed_primal_vars,
            dual_vars=unperturbed_dual_vars,
            lagrangian_gaps=lagrangian_gap_trajectory,
            gradient_norms=gradient_norms  # Return the gradient norms
        )

    def _build_vqec_result(
        self,
        energies: list[float],
        constraints: list[list[float]],
        primal_vars: np.ndarray,
        dual_vars: np.ndarray,
        lagrangian_gaps: np.ndarray,
        gradient_norms: list[float],  # Added this to return gradient norms
    ) -> VQECResult:
        """
        Builds a VQEC result.

        Args:
            energies: The energies of the optimization.
            constraints: The constraints of the optimization.
            primal_vars: The primal variables of the optimization.
            dual_vars: The dual variables of the optimization.

        Returns:
            A VQEC result.
        """
        result = VQECResult()
        result.ansatz = self._ansatz.copy()
        result.energies = np.array(energies)
        result.constraints = np.array(constraints)
        result.optimal_primal_vars = primal_vars[0]
        result.optimal_dual_vars = dual_vars
        result.vqec_iterations = len(energies) - 1
        result.lagrangian_gap_trajectory = np.array(lagrangian_gaps)
        result.gradient_norms = np.array(gradient_norms)  # Store gradient norms in result
        return result
