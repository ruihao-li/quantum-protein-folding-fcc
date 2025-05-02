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


class PerturbedPrimalDualOpt:
    """
    Class implementing the perturbed primal-dual (PPD) optimization method for the VQEC.
    """

    def __init__(
        self,
        qubit_op: SparsePauliOp,
        constr_ops: list[SparsePauliOp],
        ansatz: QuantumCircuit,
        estimator: BaseEstimator,
        gradient: BaseEstimatorGradient,
    ):
        """
        Args:
            qubit_op: The qubit operator for which the expectation value is
            minimized.
            constr_ops: The constraint operators whose the expectation values
            are constrained.
            ansatz: The quantum circuit that prepares the quantum state.
            estimator: The estimator that computes the expectation values.
            gradient: The gradient tool that computes the gradients of the
            expectation values.
        """
        self._qubit_op = qubit_op
        self._constr_ops = constr_ops
        self._ansatz = ansatz
        self._estimator = estimator
        self._gradient = gradient

    def get_expectation(
        self, circuit: QuantumCircuit, observable: SparsePauliOp, params: np.ndarray
    ) -> float:
        """
        Computes the expectation value of an observable with respect to a
        quantum state.

        Args:
            circuit: The quantum circuit that prepares the quantum state.
            observable: The observable for which the expectation value is
            computed.
            params: The parameters of the quantum circuit.

        Returns:
            The expectation value of the observable.
        """
        return self._estimator.run(circuit, observable, params).result().values[0]

    def optimize_primal_dual(
        self,
        initial_params: np.ndarray | None = None,
        initial_dual_vars: np.ndarray | None = None,
        primal_perturb_step: float = 0.05,
        dual_perturb_step: float = 0.05,
        gamma: float = 1.0,
        max_iter: int = 200,
    ) -> VQECResult:
        """
        Optimizes the primal and dual problems of the VQEC.

        Args:
            initial_params: The initial parameters of the quantum circuit.
            initial_dual_vars: The initial dual variables of the optimization.
            primal_perturb_step: The perturbation step size for the primal
            variables.
            dual_perturb_step: The perturbation step size for the dual
            variables.
            gamma: The step size parameter used for the primal-dual update.
            max_iter: The maximum number of iterations of the optimization.

        Returns:
            The result of the optimization.
        """
        # Generate the initial parameters and dual variables randomly if not provided
        initial_params = (
            initial_params
            if initial_params is not None
            else np.random.uniform(0, 2 * np.pi, self._ansatz.num_parameters)
        )
        initial_dual_vars = (
            initial_dual_vars
            if initial_dual_vars is not None
            # else np.array([0.1] * len(constr_ops))
            else np.random.uniform(0, 1, len(self._constr_ops))
        )

        unperturbed_primal_vars = np.array([initial_params])
        unperturbed_dual_vars = initial_dual_vars

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
        for i in tqdm(range(max_iter)):
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
                unperturbed_primal_vars - primal_perturb_step * primal_var_perturbations
            )
            # Perform the projection of the perturbed primal variables onto the feasible set [0, 4\pi] due to periodicity
            perturbed_primal_vars = np.mod(perturbed_primal_vars, 4 * np.pi)

            # # Orthogonal projection onto [0, 4*\pi]
            # if np.any(perturbed_primal_vars < 0):
            #     perturbed_primal_vars = np.maximum(perturbed_primal_vars, 0)
            # if np.any(perturbed_primal_vars > 4 * np.pi):
            #     perturbed_primal_vars = np.minimum(perturbed_primal_vars, 4 * np.pi)

            # Compute the perturbed dual variables
            dual_var_perturbations = np.array(constraints[-1])
            perturbed_dual_vars = np.maximum(
                unperturbed_dual_vars + dual_perturb_step * dual_var_perturbations,
                0,
            )

            # To determine the update step sizes, we follow the method proposed
            # in: M. Kallio and A. Ruszczynski, Perturbation methods for saddle
            # point computation, Tech. Rep. (International Institute for Applied
            # Systems Analysis, Laxenburg, Austria: WP-94-038, 1994).

            # Compute the Lagrangian at perturbed dual variables
            lagrangian_perturbed_dual = (
                energies[-1] + perturbed_dual_vars @ dual_var_perturbations
            )
            # Compute the Lagrangian at perturbed primal variables
            # First compute the expectation of the objective and constraint operators with the perturbed primal variables
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
                obj_exp_perturbed_primal
                + unperturbed_dual_vars @ constr_exp_perturbed_primal
            )
            # Compute the gap function
            lagrangian_gap = lagrangian_perturbed_dual - lagrangian_perturbed_primal
            lagrangian_gap_trajectory.append(lagrangian_gap)
            if np.abs(lagrangian_gap) < 1e-10:
                print(f"Gap closed after {i+1} iterations; VQEC converged.")
                break
            # d_\theta = -\nabla_\theta L(\theta, \tilde{\lambda})
            grad_primal = -(
                obj_grad + (perturbed_dual_vars * constr_grads.T).T.sum(axis=0)
            )
            # d_\lambda = \nabla_\lambda L(\tilde{\theta}, \lambda)
            grad_dual = constr_exp_perturbed_primal
            # Update step size given by: \gamma * gap / ||d_\theta||^2 + ||d_\lambda||^2
            step_size = (
                gamma
                * lagrangian_gap
                / (np.linalg.norm(grad_primal) ** 2 + np.linalg.norm(grad_dual) ** 2)
            )

            # Compute the updated primal variables
            updated_primal_vars = unperturbed_primal_vars + step_size * grad_primal
            # Projection onto the feasible set [0, 4\pi] due to periodicity
            updated_primal_vars = np.mod(updated_primal_vars, 4 * np.pi)

            # # Orthogonal projection onto [0, 4*\pi]
            # if np.any(updated_primal_vars < 0):
            #     updated_primal_vars = np.maximum(updated_primal_vars, 0)
            # if np.any(updated_primal_vars > 4 * np.pi):
            #     updated_primal_vars = np.minimum(updated_primal_vars, 4 * np.pi)

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
            primal_perturb_step=primal_perturb_step,
            dual_perturb_step=dual_perturb_step,
            primal_dual_update_step=gamma,
            energies=energies,
            constraints=constraints,
            optimal_primal_vars=unperturbed_primal_vars,
            optimal_dual_vars=unperturbed_dual_vars,
            lagrangian_gaps=lagrangian_gap_trajectory,
        )

    def _build_vqec_result(
        self,
        primal_perturb_step: float,
        dual_perturb_step: float,
        primal_dual_update_step: float,
        energies: list[float],
        constraints: list[list[float]],
        optimal_primal_vars: np.ndarray,
        optimal_dual_vars: np.ndarray,
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
        result.primal_perturb_step = primal_perturb_step
        result.dual_perturb_step = dual_perturb_step
        result.primal_dual_update_step = primal_dual_update_step
        result.energies = np.array(energies)
        result.constraints = np.array(constraints)
        result.optimal_primal_vars = optimal_primal_vars[0]
        result.optimal_dual_vars = optimal_dual_vars
        result.vqec_iterations = len(energies) - 1
        result.lagrangian_gap_trajectory = np.array(lagrangian_gaps)
        return result


class DualDecompOpt:
    """
    Class implementing the dual decomposition optimization method for the VQEC.
    """

    def __init__(
        self,
        qubit_op: SparsePauliOp,
        constr_ops: list[SparsePauliOp],
        ansatz: QuantumCircuit,
        estimator: BaseEstimator,
        optimizer: str,
    ):
        """
        Args:
            qubit_op: The qubit operator for which the expectation value is
            minimized.
            constr_ops: The constraint operators whose the expectation values
            are constrained.
            ansatz: The quantum circuit that prepares the quantum state.
            estimator: The estimator that computes the expectation values.
            optimizer: The SciPy optimizer used for the optimization.
        """
        self._qubit_op = qubit_op
        self._constr_ops = constr_ops
        self._ansatz = ansatz
        self._estimator = estimator
        self._optimizer = optimizer

    def _get_expectation(
        self, circuit: QuantumCircuit, observable: SparsePauliOp, params: np.ndarray
    ) -> float:
        """
        Computes the expectation value of an observable with respect to a
        quantum state.

        Args:
            circuit: The quantum circuit that prepares the quantum state.
            observable: The observable for which the expectation value is
            computed. params: The parameters of the quantum circuit.

        Returns:
            The expectation value of the observable.
        """
        return self._estimator.run(circuit, observable, params).result().values[0]

    def cost_function(self, primal_vars: np.ndarray, dual_vars: np.ndarray) -> float:
        """
        Computes the cost function of the VQEC.

        Args:
            primal_vars: The primal variables of the optimization.
            dual_vars: The dual variables of the optimization.

        Returns:
            The cost function value.
        """
        energy = self._get_expectation(self._ansatz, self._qubit_op, primal_vars)
        constraints = [
            self._get_expectation(self._ansatz, constr_op, primal_vars)
            for constr_op in self._constr_ops
        ]
        return energy + dual_vars @ constraints

    def optimize_primal_dual(
        self,
        initial_params: np.ndarray | None = None,
        initial_dual_vars: np.ndarray | None = None,
        dual_update_step: float = 0.1,
        max_iter: int = 200,
        opt_options: None | dict = None,
    ) -> VQECResult:
        """
        Optimizes the primal and dual problems of the VQEC.

        Args:
            initial_params: The initial parameters of the quantum circuit.
            initial_dual_vars: The initial dual variables of the optimization.
            dual_update_step: The perturbation step size for the dual variables.
            max_iter: The maximum number of iterations of the optimization.
            opt_options: The options for the SciPy optimizer

        Returns:
            The result of the optimization.
        """
        from scipy.optimize import minimize

        initial_primal_vars = (
            initial_params
            if initial_params is not None
            else np.random.uniform(0, 2 * np.pi, self._ansatz.num_parameters)
        )
        initial_dual_vars = (
            initial_dual_vars
            if initial_dual_vars is not None
            else np.array([0.1] * len(self._constr_ops))
        )
        self.energy_trajectory = [
            self._get_expectation(self._ansatz, self._qubit_op, initial_primal_vars)
        ]
        self.constraints_trajectory = [
            [
                self._get_expectation(self._ansatz, constr_op, initial_primal_vars)
                for constr_op in self._constr_ops
            ]
        ]

        self.primal_vars = initial_primal_vars
        self.dual_vars = initial_dual_vars
        self.dual_vars_history = [self.dual_vars]

        if opt_options is None:
            opt_options = {}

        for i in tqdm(range(max_iter)):
            # Full minimization of the Lagrangian w.r.t. primal variables
            primal_result = minimize(
                fun=self.cost_function,
                x0=self.primal_vars,
                args=(self.dual_vars),
                method=self._optimizer,
                options=opt_options,
            )
            self.primal_vars = primal_result.x

            # Update the dual variables
            dual_update_step_i = dual_update_step / (i + 1)
            constraints = np.array(
                [
                    self._get_expectation(self._ansatz, constr_op, self.primal_vars)
                    for constr_op in self._constr_ops
                ]
            )
            updated_dual_vars = np.maximum(
                self.dual_vars + dual_update_step_i * constraints, 0
            )
            self.dual_vars = updated_dual_vars

            self.energy_trajectory.append(
                self._get_expectation(self._ansatz, self._qubit_op, self.primal_vars)
            )
            self.constraints_trajectory.append(constraints)
            self.dual_vars_history.append(updated_dual_vars)

            # Check for convergence
            tol = 1e-5
            if (
                np.abs(self.energy_trajectory[-1] - self.energy_trajectory[-2])
                / np.abs(self.energy_trajectory[-2])
                < tol
            ):
                print(f"Converged after {i+1} iterations")
                break

        return self._build_vqec_result()

    def _build_vqec_result(self) -> VQECResult:
        """
        Builds a VQEC result.

        Returns:
            A VQEC result.
        """
        result = VQECResult()
        result.ansatz = self._ansatz.copy()
        result.energies = np.array(self.energy_trajectory)
        result.constraints = np.array(self.constraints_trajectory)
        result.dual_vars_history = np.array(self.dual_vars_history)
        result.optimal_primal_vars = self.primal_vars
        result.optimal_dual_vars = self.dual_vars
        result.vqec_iterations = len(self.energy_trajectory) - 1
        return result
