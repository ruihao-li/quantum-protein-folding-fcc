"""Variational Quantum Eigensolver with Constraints (VQEC) module. Based on
arXiv:2311.08502. This version contains an automatic selection of the update
step size based on the method proposed in: M. Kallio and A. Ruszczynski,
Perturbation methods for saddle point computation, Tech. Rep. (International
Institute for Applied Systems Analysis, Laxenburg, Austria: WP-94-038, 1994)."""

from __future__ import annotations
from .vqec import VQECResult
from qiskit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp
from qiskit.primitives import BaseEstimator
from qiskit_algorithms.gradients import BaseEstimatorGradient
import numpy as np
from tqdm import tqdm


class VQECOpt:

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
            initial_params
            if initial_params is not None
            else np.random.uniform(0, 2 * np.pi, ansatz.num_parameters)
        )
        self._initial_dual_vars = (
            initial_dual_vars
            if initial_dual_vars is not None
            else np.array([0.1] * len(constr_ops))
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
        for i in tqdm(range(self._max_iter)):
            # Compute gradients for the objective and constraint operators
            obj_grad = np.array(
                (
                    self._gradient.run(
                        circuits=self._ansatz,
                        observables=self._qubit_op,
                        parameter_values=unperturbed_primal_vars,
                    )
                    .result()
                    .gradients
                )
            )
            constr_grads = np.array(
                [
                    self._gradient.run(
                        circuits=self._ansatz,
                        observables=constr_op,
                        parameter_values=unperturbed_primal_vars,
                    )
                    .result()
                    .gradients
                    for constr_op in self._constr_ops
                ]
            )

            # Compute the perturbed primal variables
            # unperturbed_dual_vars dim: (n_constraints,)
            # constr_grads dim: (n_constraints, 1, n_parameters)
            # obj_grad dim: (1, n_parameters)
            primal_var_perturbations = (unperturbed_dual_vars * constr_grads.T).T.sum(
                axis=0
            ) + obj_grad
            perturbed_primal_vars = (
                unperturbed_primal_vars
                - self._primal_perturb_step * primal_var_perturbations
            )

            # Compute the perturbed dual variables
            dual_var_perturbations = np.array(
                [
                    self.get_expectation(
                        self._ansatz, constr_op, unperturbed_primal_vars
                    )
                    for constr_op in self._constr_ops
                ]
            )
            perturbed_dual_vars = np.maximum(
                unperturbed_dual_vars
                + self._dual_perturb_step * dual_var_perturbations,
                0,
            )

            # To determine the update step sizes, we follow the method proposed
            # in: M. Kallio and A. Ruszczynski, Perturbation methods for saddle
            # point computation, Tech. Rep. (International Institute for Applied
            # Systems Analysis, Laxenburg, Austria: WP-94-038, 1994).

            # Compute the Lagrangian at perturbed dual variables
            lagrangian_perturbed_dual = self.get_expectation(
                self._ansatz, self._qubit_op, unperturbed_primal_vars
            ) + perturbed_dual_vars @ [
                self.get_expectation(self._ansatz, constr_op, unperturbed_primal_vars)
                for constr_op in self._constr_ops
            ]
            # Compute the Lagrangian at perturbed primal variables
            lagrangian_perturbed_primal = self.get_expectation(
                self._ansatz, self._qubit_op, perturbed_primal_vars
            ) + unperturbed_dual_vars @ [
                self.get_expectation(self._ansatz, constr_op, perturbed_primal_vars)
                for constr_op in self._constr_ops
            ]
            lagrangian_gap = lagrangian_perturbed_dual - lagrangian_perturbed_primal
            lagrangian_gap_trajectory.append(lagrangian_gap)
            if lagrangian_gap < 1e-10:
                print(f"Gap closed after {i+1} iterations; VQEC converged.")
                break
            # d_\theta = -\nabla_\theta L(\theta, \tilde{\lambda})
            grad_primal = -(
                obj_grad + (perturbed_dual_vars * constr_grads.T).T.sum(axis=0)
            )
            # d_\lambda = \nabla_\lambda L(\tilde{\theta}, \lambda)
            grad_dual = np.array(
                [
                    self.get_expectation(self._ansatz, constr_op, perturbed_primal_vars)
                    for constr_op in self._constr_ops
                ]
            )
            # Update step size given by: \gamma * gap / ||d_\theta||^2 + ||d_\lambda||^2
            step_size = (
                self._gamma
                * lagrangian_gap
                / (np.linalg.norm(grad_primal) ** 2 + np.linalg.norm(grad_dual) ** 2)
            )

            # Compute the updated primal variables
            updated_primal_vars = unperturbed_primal_vars + step_size * grad_primal

            # Compute the updated dual variables
            updated_dual_vars = np.maximum(
                unperturbed_dual_vars + step_size * grad_dual, 0
            )

            # Compute expectation values with updated primal variables
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
        )

    def _build_vqec_result(
        self,
        energies: list[float],
        constraints: list[list[float]],
        primal_vars: np.ndarray,
        dual_vars: np.ndarray,
        lagrangian_gaps: np.ndarray,
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
        return result


# class VQECResult:
#     def __init__(
#         self,
#         ansatz: QuantumCircuit | None = None,
#         energies: np.ndarray | None = None,
#         constraints: np.ndarray | None = None,
#         optimal_primal_vars: np.ndarray | None = None,
#         optimal_dual_vars: np.ndarray | None = None,
#         vqec_iterations: int | None = None,
#     ):
#         """
#         Args:
#             ansatz: The quantum circuit that prepares the quantum state.
#             energies: The energies of the optimization.
#             constraints: The constraints of the optimization.
#             optimal_primal_vars: The optimal primal variables of the optimization.
#             optimal_dual_vars: The optimal dual variables of the optimization.
#             vqec_iterations: The number of iterations of the optimization.
#         """
#         self._ansatz = ansatz
#         self._energies = energies
#         self._constraints = constraints
#         self._optimal_primal_vars = optimal_primal_vars
#         self._optimal_dual_vars = optimal_dual_vars
#         self._vqec_iterations = vqec_iterations

#     @property
#     def ansatz(self) -> QuantumCircuit:
#         """Returns the quantum circuit that prepares the quantum state."""
#         return self._ansatz

#     @ansatz.setter
#     def ansatz(self, ansatz: QuantumCircuit):
#         """Sets the quantum circuit that prepares the quantum state."""
#         self._ansatz = ansatz

#     @property
#     def energies(self) -> np.ndarray:
#         """Returns the energies of the optimization."""
#         return self._energies

#     @energies.setter
#     def energies(self, energies: np.ndarray):
#         """Sets the energies of the optimization."""
#         self._energies = energies

#     @property
#     def constraints(self) -> np.ndarray:
#         """Returns the constraints of the optimization."""
#         return self._constraints

#     @constraints.setter
#     def constraints(self, constraints: np.ndarray):
#         """Sets the constraints of the optimization. Note that the shape of the list is (number of iterations, number of constraints)."""
#         self._constraints = constraints

#     @property
#     def optimal_primal_vars(self) -> np.ndarray:
#         """Returns the optimal primal variables of the optimization."""
#         return self._optimal_primal_vars

#     @optimal_primal_vars.setter
#     def optimal_primal_vars(self, optimal_primal_vars: np.ndarray):
#         """Sets the optimal primal variables of the optimization."""
#         self._optimal_primal_vars = optimal_primal_vars

#     @property
#     def optimal_dual_vars(self) -> np.ndarray:
#         """Returns the optimal dual variables of the optimization."""
#         return self._optimal_dual_vars

#     @optimal_dual_vars.setter
#     def optimal_dual_vars(self, optimal_dual_vars: np.ndarray):
#         """Sets the optimal dual variables of the optimization."""
#         self._optimal_dual_vars = optimal_dual_vars

#     @property
#     def vqec_iterations(self) -> int:
#         """Returns the number of iterations of the optimization."""
#         return self._vqec_iterations

#     @vqec_iterations.setter
#     def vqec_iterations(self, vqec_iterations: int):
#         """Sets the number of iterations of the optimization."""
#         self._vqec_iterations = vqec_iterations
