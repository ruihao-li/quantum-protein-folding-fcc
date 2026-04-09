# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

"""
Variational Quantum Eigensolver with Constraints (VQEC) module. Based on
arXiv:2311.08502.
This version contains an option for automatic selection of the update step size
based on the method proposed in: M. Kallio and A. Ruszczynski, Perturbation
methods for saddle point computation, Tech. Rep. (International Institute for
Applied Systems Analysis, Laxenburg, Austria: WP-94-038, 1994).
"""

from __future__ import annotations
from qiskit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp
from qiskit.primitives import BaseEstimatorV2
from qiskit_algorithms.gradients import BaseEstimatorGradient
import numpy as np
from tqdm import tqdm


def _normalize_primal_params(params: np.ndarray, num_parameters: int) -> np.ndarray:
    """Normalize a single variational parameter set to shape ``(1, num_parameters)``."""
    values = np.asarray(params, dtype=float)
    if values.ndim == 0:
        raise ValueError("initial_params must contain at least one parameter value.")
    if values.ndim == 1:
        values = values.reshape(1, -1)
    elif values.ndim != 2:
        raise ValueError(
            "initial_params must be one parameter vector or one batched parameter row."
        )
    if values.shape[0] != 1:
        raise ValueError("Exactly one parameter row is supported for VQEC optimization.")
    if values.shape[1] != num_parameters:
        raise ValueError(
            f"Expected {num_parameters} ansatz parameters, got {values.shape[1]}."
        )
    return values.copy()


def _normalize_dual_vars(dual_vars: np.ndarray, num_constraints: int) -> np.ndarray:
    """Normalize dual variables to a flat vector of length ``num_constraints``."""
    values = np.asarray(dual_vars, dtype=float).reshape(-1)
    if values.size != num_constraints:
        raise ValueError(
            f"Expected {num_constraints} dual variables, got {values.size}."
        )
    return values.copy()


def _first_expectation_value(result) -> float:
    """Return the first expectation value from an Estimator V2 PUB result."""
    evs = np.asarray(result[0].data.evs)
    if evs.size == 0:
        raise ValueError("Estimator returned no expectation values.")
    return float(evs.reshape(-1)[0])


def _validate_operator_shapes(
    qubit_op: SparsePauliOp,
    constr_ops: list[SparsePauliOp],
    ansatz: QuantumCircuit,
) -> None:
    """Validate that the ansatz and all observables act on the same qubits."""
    if qubit_op.num_qubits != ansatz.num_qubits:
        raise ValueError(
            "The ansatz and objective operator must act on the same number of qubits: "
            f"{ansatz.num_qubits} != {qubit_op.num_qubits}."
        )
    for idx, constr_op in enumerate(constr_ops):
        if constr_op.num_qubits != ansatz.num_qubits:
            raise ValueError(
                "Constraint operator qubit counts must match the ansatz qubit count. "
                f"Constraint {idx} uses {constr_op.num_qubits} qubits while the ansatz uses {ansatz.num_qubits}."
            )


class PerturbedPrimalDualOpt:
    """
    Class implementing the perturbed primal-dual (PPD) optimization method for
    the VQEC. Based on Le & Kekatos, arXiv:2311.08502.
    """

    def __init__(
        self,
        qubit_op: SparsePauliOp,
        constr_ops: list[SparsePauliOp],
        ansatz: QuantumCircuit,
        estimator: BaseEstimatorV2,
        gradient: BaseEstimatorGradient | None = None,
    ):
        """
        Args:
            qubit_op: The qubit operator for which the expectation value is
            minimized.
            constr_ops: The constraint operators whose expectation values
            are constrained.
            ansatz: The quantum circuit that prepares the quantum state.
            estimator: The estimator that computes the expectation values.
            gradient: Optional gradient tool that computes the gradients of the
            expectation values. When the provided gradient helper is
            incompatible with Estimator V2, the optimizer falls back to a
            finite-difference gradient computed through ``estimator``.
        """
        _validate_operator_shapes(qubit_op, constr_ops, ansatz)
        self._qubit_op = qubit_op
        self._constr_ops = constr_ops
        self._ansatz = ansatz
        self._estimator = estimator
        self._gradient = gradient
        self._use_finite_differences = gradient is None

    def _get_gradient(
        self, observable: SparsePauliOp, params: np.ndarray, eps: float = 1e-6
    ) -> np.ndarray:
        """Return a gradient row compatible with the current estimator backend."""
        parameter_values = _normalize_primal_params(params, self._ansatz.num_parameters)
        if self._ansatz.num_parameters == 0:
            return np.zeros((1, 0), dtype=float)

        if not self._use_finite_differences and self._gradient is not None:
            try:
                gradients = (
                    self._gradient.run(
                        circuits=self._ansatz,
                        observables=observable,
                        parameter_values=parameter_values,
                    )
                    .result()
                    .gradients
                )
                return np.asarray(gradients, dtype=float).reshape(1, -1)
            except TypeError as exc:
                if "EstimatorV2.run" not in str(exc):
                    raise
                self._use_finite_differences = True

        base_params = parameter_values.reshape(-1)
        gradient = np.zeros((1, self._ansatz.num_parameters), dtype=float)
        for param_idx in range(self._ansatz.num_parameters):
            shifted_plus = base_params.copy()
            shifted_minus = base_params.copy()
            shifted_plus[param_idx] += eps
            shifted_minus[param_idx] -= eps
            gradient[0, param_idx] = (
                self.get_expectation(self._ansatz, observable, shifted_plus)
                - self.get_expectation(self._ansatz, observable, shifted_minus)
            ) / (2 * eps)
        return gradient

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
        parameter_values = _normalize_primal_params(params, circuit.num_parameters)
        result = self._estimator.run([(circuit, observable, parameter_values)]).result()
        return _first_expectation_value(result)

    def optimize_primal_dual(
        self,
        initial_params: np.ndarray | None = None,
        initial_dual_vars: np.ndarray | None = None,
        primal_perturb_step: float = 0.05,
        dual_perturb_step: float = 0.05,
        gamma: float = 1.0,
        auto_update_step: bool = True,
        max_iter: int = 200,
    ) -> dict:
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
            auto_update_step: Whether to automatically update the step size
            based on the method proposed in: M. Kallio and A. Ruszczynski,
            Perturbation methods for saddle point computation, Tech. Rep.
            (International Institute for Applied Systems Analysis, Laxenburg,
            Austria: WP-94-038, 1994). If False, the step size is dynamically
            adjusted based on an exponential decay of gamma.
            max_iter: The maximum number of iterations of the optimization.

        Returns:
            The result of the optimization.
        """
        # Generate the initial parameters and dual variables randomly if not provided
        initial_params = (
            _normalize_primal_params(initial_params, self._ansatz.num_parameters)
            if initial_params is not None
            else np.random.uniform(0, 2 * np.pi, (1, self._ansatz.num_parameters))
        )
        initial_dual_vars = (
            _normalize_dual_vars(initial_dual_vars, len(self._constr_ops))
            if initial_dual_vars is not None
            # else np.array([0.1] * len(constr_ops))
            else np.random.uniform(0, 1, len(self._constr_ops))
        )

        unperturbed_primal_vars = initial_params.copy()
        unperturbed_dual_vars = initial_dual_vars.copy()

        energies = [
            self.get_expectation(self._ansatz, self._qubit_op, unperturbed_primal_vars)
        ]
        constraints = [
            [
                self.get_expectation(self._ansatz, constr_op, unperturbed_primal_vars)
                for constr_op in self._constr_ops
            ]
        ]
        update_step_history = []
        for i in tqdm(range(max_iter)):
            # Compute gradients for the objective and constraint operators
            obj_grad = self._get_gradient(self._qubit_op, unperturbed_primal_vars)
            constr_grads = np.array(
                [
                    self._get_gradient(constr_op, unperturbed_primal_vars)
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
            # Perform the projection of the perturbed primal variables onto the feasible set [0, 2\pi] due to periodicity
            perturbed_primal_vars = np.mod(perturbed_primal_vars, 2 * np.pi)

            # # Orthogonal projection onto [0, 2*\pi]
            # if np.any(perturbed_primal_vars < 0):
            #     perturbed_primal_vars = np.maximum(perturbed_primal_vars, 0)
            # if np.any(perturbed_primal_vars > 2 * np.pi):
            #     perturbed_primal_vars = np.minimum(perturbed_primal_vars, 2 * np.pi)

            # Compute the perturbed dual variables
            dual_var_perturbations = np.array(constraints[-1])
            perturbed_dual_vars = np.maximum(
                unperturbed_dual_vars + dual_perturb_step * dual_var_perturbations,
                0,
            )

            if auto_update_step:
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
                        self.get_expectation(
                            self._ansatz, constr_op, perturbed_primal_vars
                        )
                        for constr_op in self._constr_ops
                    ]
                )
                lagrangian_perturbed_primal = (
                    obj_exp_perturbed_primal
                    + unperturbed_dual_vars @ constr_exp_perturbed_primal
                )
                # Compute the gap function
                lagrangian_gap = lagrangian_perturbed_dual - lagrangian_perturbed_primal
                # if np.abs(lagrangian_gap) < 1e-10:
                #     print(f"Gap closed after {i+1} iterations; VQEC converged.")
                #     break
                # d_\theta = -\nabla_\theta L(\theta, \tilde{\lambda})
                grad_primal = -(
                    obj_grad + (perturbed_dual_vars * constr_grads.T).T.sum(axis=0)
                )
                # d_\lambda = \nabla_\lambda L(\tilde{\theta}, \lambda)
                grad_dual = constr_exp_perturbed_primal
                # Update step size given by: \gamma * gap / (||d_\theta||^2 + ||d_\lambda||^2)
                grad_norm_sq = (
                    np.linalg.norm(grad_primal) ** 2
                    + np.linalg.norm(grad_dual) ** 2
                )
                if grad_norm_sq <= 1e-12:
                    step_size = 0.0
                else:
                    step_size = float(gamma * lagrangian_gap / grad_norm_sq)
            else:
                # d_\theta = -\nabla_\theta L(\theta, \tilde{\lambda})
                grad_primal = -(
                    obj_grad + (perturbed_dual_vars * constr_grads.T).T.sum(axis=0)
                )
                # d_\lambda = \nabla_\lambda L(\tilde{\theta}, \lambda)
                constr_exp_perturbed_primal = np.array(
                    [
                        self.get_expectation(
                            self._ansatz, constr_op, perturbed_primal_vars
                        )
                        for constr_op in self._constr_ops
                    ]
                )
                grad_dual = constr_exp_perturbed_primal
                step_size = gamma * 0.99**i
            # Store the update step size
            update_step_history.append(step_size)

            # Compute the updated primal variables
            updated_primal_vars = unperturbed_primal_vars + step_size * grad_primal
            # Projection onto the feasible set [0, 2\pi] due to periodicity
            updated_primal_vars = np.mod(updated_primal_vars, 2 * np.pi)

            # # Orthogonal projection onto [0, 2*\pi]
            # if np.any(updated_primal_vars < 0):
            #     updated_primal_vars = np.maximum(updated_primal_vars, 0)
            # if np.any(updated_primal_vars > 2 * np.pi):
            #     updated_primal_vars = np.minimum(updated_primal_vars, 2 * np.pi)

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

            # Store the latest iterate before the convergence check so the
            # returned optimum matches the final evaluated point.
            unperturbed_primal_vars = updated_primal_vars
            unperturbed_dual_vars = updated_dual_vars

            # Check for convergence
            tol = 1e-5
            denominator = max(np.abs(energies[-2]), 1e-12)
            if np.abs(energies[-1] - energies[-2]) / denominator < tol:
                print(f"Converged after {i+1} iterations")
                break

        result = {
            "primal_perturb_step": primal_perturb_step,
            "dual_perturb_step": dual_perturb_step,
            "update_step_history": update_step_history,
            "energy_history": energies,
            "constraints_history": constraints,
            "optimal_primal_vars": np.asarray(unperturbed_primal_vars).reshape(-1).copy(),
            "optimal_dual_vars": np.asarray(unperturbed_dual_vars).reshape(-1).copy(),
        }

        return result


class OptimisticGDAOpt:
    """
    Class implementing the optimistic gradient descent ascent (OGDA)
    optimization method for the VQEC. Based on Mokhtari et al.,
    arXiv:1901.08511.
    """

    def __init__(
        self,
        qubit_op: SparsePauliOp,
        constr_ops: list[SparsePauliOp],
        ansatz: QuantumCircuit,
        estimator: BaseEstimatorV2,
        gradient: BaseEstimatorGradient | None = None,
    ):
        """
        Args:
            qubit_op: The qubit operator for which the expectation value is
            minimized.
            constr_ops: The constraint operators whose expectation values
            are constrained.
            ansatz: The quantum circuit that prepares the quantum state.
            estimator: The estimator that computes the expectation values.
            gradient: Optional gradient tool that computes the gradients of the
            expectation values. When the provided gradient helper is
            incompatible with Estimator V2, the optimizer falls back to a
            finite-difference gradient computed through ``estimator``.
        """
        _validate_operator_shapes(qubit_op, constr_ops, ansatz)
        self._qubit_op = qubit_op
        self._constr_ops = constr_ops
        self._ansatz = ansatz
        self._estimator = estimator
        self._gradient = gradient
        self._use_finite_differences = gradient is None

    def _get_gradient(
        self, observable: SparsePauliOp, params: np.ndarray, eps: float = 1e-6
    ) -> np.ndarray:
        """Return a gradient row compatible with the current estimator backend."""
        parameter_values = _normalize_primal_params(params, self._ansatz.num_parameters)
        if self._ansatz.num_parameters == 0:
            return np.zeros((1, 0), dtype=float)

        if not self._use_finite_differences and self._gradient is not None:
            try:
                gradients = (
                    self._gradient.run(
                        circuits=self._ansatz,
                        observables=observable,
                        parameter_values=parameter_values,
                    )
                    .result()
                    .gradients
                )
                return np.asarray(gradients, dtype=float).reshape(1, -1)
            except TypeError as exc:
                if "EstimatorV2.run" not in str(exc):
                    raise
                self._use_finite_differences = True

        base_params = parameter_values.reshape(-1)
        gradient = np.zeros((1, self._ansatz.num_parameters), dtype=float)
        for param_idx in range(self._ansatz.num_parameters):
            shifted_plus = base_params.copy()
            shifted_minus = base_params.copy()
            shifted_plus[param_idx] += eps
            shifted_minus[param_idx] -= eps
            gradient[0, param_idx] = (
                self.get_expectation(self._ansatz, observable, shifted_plus)
                - self.get_expectation(self._ansatz, observable, shifted_minus)
            ) / (2 * eps)
        return gradient

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
        parameter_values = _normalize_primal_params(params, circuit.num_parameters)
        result = self._estimator.run([(circuit, observable, parameter_values)]).result()
        return _first_expectation_value(result)

    def optimize_primal_dual(
        self,
        initial_params: np.ndarray | None = None,
        initial_dual_vars: np.ndarray | None = None,
        learning_rate: float = 0.01,
        schedule: str = "constant",
        alpha: float = 1.0,
        beta: float = 1.0,
        max_iter: int = 200,
        lr_options: None | dict = None,
    ) -> dict:
        r"""
        Optimizes the primal and dual parameters of the VQEC, based on the following updates:
        $$
        \theta_{k+1} = \theta_k - (\alpha + \beta)\eta_k \nabla_\theta L(\theta_k, \lambda_k) + \beta \eta_k \nabla_\theta L(\theta_{k-1}, \lambda_{k-1}) \\
        \lambda_{k+1} = \lambda_k + (\alpha + \beta)\eta_k \nabla_\lambda L(\theta_k, \lambda_k) - \beta \eta_k \nabla_\lambda L(\theta_{k-1}, \lambda_{k-1})
        $$

        Args:
            initial_params: The initial parameters of the quantum circuit.
            initial_dual_vars: The initial dual variables of the optimization.
            learning_rate: The learning rate for the optimization.
            schedule: The schedule for the learning rate. Options: "constant",
            "inverse_linear", "exponential".
            alpha: The parameter for the optimistic update.
            beta: The parameter for the negative momentum.
            max_iter: The maximum number of iterations of the optimization.
            lr_options: The options for the learning rate schedule. Options:
            "decay_rate". Default is None.

        Returns:
            The result of the optimization as a dictionary.
        """
        # Generate the initial parameters and dual variables randomly if not provided
        initial_params = (
            _normalize_primal_params(initial_params, self._ansatz.num_parameters)
            if initial_params is not None
            else np.random.uniform(0, 2 * np.pi, (1, self._ansatz.num_parameters))
        )
        initial_dual_vars = (
            _normalize_dual_vars(initial_dual_vars, len(self._constr_ops))
            if initial_dual_vars is not None
            # else np.array([0.1] * len(constr_ops))
            else np.random.uniform(0, 1, len(self._constr_ops))
        )

        primal_vars_list = [initial_params.copy()]
        dual_vars_list = [initial_dual_vars.copy()]

        energies = [
            self.get_expectation(self._ansatz, self._qubit_op, primal_vars_list[0])
        ]
        constraints = [
            [
                self.get_expectation(self._ansatz, constr_op, primal_vars_list[0])
                for constr_op in self._constr_ops
            ]
        ]
        obj_grads_list = []
        constr_grads_list = []
        dynamic_lr = learning_rate
        for i in tqdm(range(max_iter)):
            # Compute gradients for the objective and constraint operators with current primal variables
            obj_grad = self._get_gradient(self._qubit_op, primal_vars_list[-1])
            constr_grads = np.array(
                [
                    self._get_gradient(constr_op, primal_vars_list[-1])
                    for constr_op in self._constr_ops
                ]
            )
            obj_grads_list.append(obj_grad)
            constr_grads_list.append(constr_grads)

            # Compute the primal update at current step
            # dual_vars dim: (n_constraints,)
            # constr_grads dim: (n_constraints, 1, n_parameters)
            # obj_grad dim: (1, n_parameters)
            primal_var_update_current = (
                dual_vars_list[-1] * constr_grads_list[-1].T
            ).T.sum(axis=0) + obj_grads_list[-1]
            # Compute the primal update at previous step if i > 0 (for the memory-based negative momentum)
            if i > 0:
                primal_var_update_prev = (
                    dual_vars_list[-2] * constr_grads_list[-2].T
                ).T.sum(axis=0) + obj_grads_list[-2]
            else:
                primal_var_update_prev = np.zeros_like(primal_var_update_current)

            updated_primal_vars = (
                primal_vars_list[-1]
                - (alpha + beta) * dynamic_lr * primal_var_update_current
                + beta * dynamic_lr * primal_var_update_prev
            )
            # Perform the projection of the perturbed primal variables onto the feasible set [0, 2*\pi] due to periodicity
            updated_primal_vars = np.mod(updated_primal_vars, 2 * np.pi)
            primal_vars_list.append(updated_primal_vars)

            # # Orthogonal projection onto [0, 2*\pi]
            # if np.any(updated_primal_vars < 0):
            #     updated_primal_vars = np.maximum(updated_primal_vars, 0)
            # if np.any(updated_primal_vars > 2 * np.pi):
            #     updated_primal_vars = np.minimum(updated_primal_vars, 2 * np.pi)

            # Compute the dual update at current step
            dual_var_update_current = np.array(constraints[-1])
            # Compute the dual update at previous step if i > 0 (for the memory-based negative momentum)
            if i > 0:
                dual_var_update_prev = np.array(constraints[-2])
            else:
                dual_var_update_prev = np.zeros_like(dual_var_update_current)
            updated_dual_vars = (
                dual_vars_list[-1]
                + (alpha + beta) * dynamic_lr * dual_var_update_current
                - beta * dynamic_lr * dual_var_update_prev
            )
            # Make sure the dual variables are non-negative
            updated_dual_vars = np.maximum(
                updated_dual_vars,
                0,
            )
            dual_vars_list.append(updated_dual_vars)

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
            tol = 1e-7
            denominator = max(np.abs(energies[-2]), 1e-12)
            if np.abs(energies[-1] - energies[-2]) / denominator < tol:
                print(f"Converged after {i+1} iterations")
                break

            # Update the step size based on the learning rate
            if schedule == "constant":
                dynamic_lr = learning_rate
            elif schedule == "inverse_linear":
                if lr_options is not None and "decay_rate" in lr_options:
                    decay_rate = lr_options["decay_rate"]
                    dynamic_lr = learning_rate * 1 / (1 + decay_rate * i)
                else:
                    dynamic_lr = learning_rate * 1 / (1 + i / max_iter)
            elif schedule == "exponential":
                if lr_options is not None and "decay_rate" in lr_options:
                    decay_rate = lr_options["decay_rate"]
                    dynamic_lr = learning_rate * np.exp(-decay_rate * i)
                else:
                    dynamic_lr = learning_rate * np.exp(-i / max_iter)
            else:
                raise ValueError(
                    "Invalid schedule type. Choose 'constant', 'inverse_linear', or 'exponential'."
                )

        result = {
            "learning_rate": learning_rate,
            "schedule": schedule,
            "alpha": alpha,
            "beta": beta,
            "energy_history": energies,
            "constraints_history": constraints,
            "optimal_primal_vars": np.asarray(primal_vars_list[-1]).reshape(-1).copy(),
            "optimal_dual_vars": np.asarray(dual_vars_list[-1]).reshape(-1).copy(),
        }

        return result
