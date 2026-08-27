"""Expectation, gradient, and PDP tests for chance-constrained VQEC."""

from __future__ import annotations

import unittest

import numpy as np
from qiskit.circuit.library import real_amplitudes
from qiskit.quantum_info import Statevector
from qiskit_aer.primitives import SamplerV2 as Sampler

import fcc
import vqe
from fcc import build_turn_only_fcc_model
from vqe import ChanceConstrainedVQEC


class UnusedSampler:
    """Satisfy constructor validation for exact Statevector tests."""

    def run(self, *args, **kwargs):  # pragma: no cover - must remain unused
        raise AssertionError("The exact test adapter must not sample")


def exact_primitives(model, ansatz, parameters: np.ndarray) -> np.ndarray:
    """Evaluate model primitives from exact Statevector probabilities."""

    batch = np.atleast_2d(np.asarray(parameters, dtype=float))
    rows = []
    for values in batch:
        state = Statevector(ansatz.assign_parameters(values))
        rows.append(model.primitive_values(state.probabilities_dict())[0])
    return np.asarray(rows)


class ExactChanceConstrainedVQEC(ChanceConstrainedVQEC):
    """Use exact probabilities while exercising production shift/PDP hooks."""

    def primitive_values(self, parameters: np.ndarray) -> np.ndarray:
        return exact_primitives(self.problem, self._ansatz, parameters)


class ChanceConstrainedVQECTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = build_turn_only_fcc_model("ACDE", objective_scale=20.0)
        cls.ansatz = real_amplitudes(
            cls.model.num_qubits, reps=1, entanglement="linear"
        )
        cls.parameters = np.linspace(
            0.13, 1.37, cls.ansatz.num_parameters, dtype=float
        )

    def test_public_package_exports_revised_and_legacy_apis(self) -> None:
        for name in (
            "ProteinFoldingProblem",
            "ProteinSolver",
            "TurnOnlyFCCModel",
            "build_turn_only_fcc_model",
        ):
            self.assertIn(name, fcc.__all__)
        for name in (
            "PerturbedPrimalDualOpt",
            "OptimisticGDAOpt",
            "ChanceConstrainedVQEC",
        ):
            self.assertIn(name, vqe.__all__)

    def test_constraint_limits_are_explicit_and_validated(self) -> None:
        with self.assertRaises(ValueError):
            ChanceConstrainedVQEC(
                self.model,
                self.ansatz,
                UnusedSampler(),
                constraint_limits=np.empty(0),
            )
        with self.assertRaises(ValueError):
            ChanceConstrainedVQEC(
                self.model,
                self.ansatz,
                UnusedSampler(),
                constraint_limits=np.asarray([1.1]),
            )

        measured = self.ansatz.copy()
        measured.measure_all()
        with self.assertRaises(ValueError):
            ChanceConstrainedVQEC(
                self.model,
                measured,
                UnusedSampler(),
                constraint_limits=np.asarray([0.0]),
            )

    def test_sampler_primitives_agree_with_exact_probabilities(self) -> None:
        shots = 40_000
        solver = ChanceConstrainedVQEC(
            self.model,
            self.ansatz,
            Sampler(default_shots=shots, seed=7103),
            constraint_limits=np.asarray([0.05]),
            shots=shots,
        )
        sampled = solver.primitive_values(self.parameters)[0]
        exact = exact_primitives(self.model, self.ansatz, self.parameters)[0]
        np.testing.assert_allclose(sampled, exact, atol=0.025, rtol=0.0)

    def test_parameter_shift_matches_finite_difference_exactly(self) -> None:
        solver = ExactChanceConstrainedVQEC(
            self.model,
            self.ansatz,
            UnusedSampler(),
            constraint_limits=np.asarray([0.07]),
        )
        objective_gradient, constraint_gradients = solver._evaluate_gradients(
            self.parameters
        )
        parameter_shift = np.vstack(
            (objective_gradient, constraint_gradients[:, 0, :])
        )

        step = 1e-6
        finite_difference = np.empty_like(parameter_shift)
        for index in range(self.ansatz.num_parameters):
            plus = self.parameters.copy()
            minus = self.parameters.copy()
            plus[index] += step
            minus[index] -= step
            plus_values = exact_primitives(self.model, self.ansatz, plus)[0]
            minus_values = exact_primitives(self.model, self.ansatz, minus)[0]
            finite_difference[:, index] = (plus_values - minus_values) / (2 * step)

        np.testing.assert_allclose(
            parameter_shift, finite_difference, atol=2e-7, rtol=0.0
        )

    def test_chance_adapter_reuses_pdp_result_contract(self) -> None:
        delta = np.asarray([0.07])
        solver = ExactChanceConstrainedVQEC(
            self.model,
            self.ansatz,
            UnusedSampler(),
            constraint_limits=delta,
        )
        exact = exact_primitives(self.model, self.ansatz, self.parameters)[0]
        objective, residuals = solver._evaluate_primitives(self.parameters)
        self.assertAlmostEqual(objective, exact[0])
        np.testing.assert_allclose(residuals, exact[1:] - delta)

        result = solver.optimize_primal_dual(
            initial_params=self.parameters,
            initial_dual_vars=np.zeros(self.model.constraint_count),
            primal_perturb_step=0.01,
            dual_perturb_step=0.01,
            gamma=0.01,
            auto_update_step=False,
            max_iter=1,
        )
        self.assertEqual(
            set(result),
            {
                "primal_perturb_step",
                "dual_perturb_step",
                "update_step_history",
                "energy_history",
                "constraints_history",
                "optimal_primal_vars",
                "optimal_dual_vars",
            },
        )
        np.testing.assert_allclose(result["constraints_history"][0], residuals)
        self.assertTrue(np.all(result["optimal_primal_vars"] >= 0.0))
        self.assertTrue(np.all(result["optimal_primal_vars"] < 2.0 * np.pi))
        self.assertTrue(np.all(result["optimal_dual_vars"] >= 0.0))


if __name__ == "__main__":
    unittest.main()
