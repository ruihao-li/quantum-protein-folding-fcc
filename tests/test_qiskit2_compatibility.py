"""Compatibility checks for the public Qiskit 2.x API surface."""

from __future__ import annotations

import unittest

import numpy as np
from qiskit.circuit.library import real_amplitudes
from qiskit.primitives import StatevectorEstimator
from qiskit.quantum_info import SparsePauliOp
from qiskit_algorithms.gradients import ParamShiftEstimatorGradient

from fcc.fcc_protein_folding_problem import ProteinFoldingProblem
from fcc.utils import compose_IZ_ops
from vqe import ChanceConstrainedVQEC, OptimisticGDAOpt, PerturbedPrimalDualOpt


class _MeasurementData:
    def get_counts(self) -> dict[str, int]:
        return {"001": 2, "111": 5}


class _PubResult:
    class Data:
        meas = _MeasurementData()

    data = Data()


class _PrimitiveResult:
    def __getitem__(self, index: int) -> _PubResult:
        if index != 0:
            raise IndexError(index)
        return _PubResult()


class Qiskit2CompatibilityTests(unittest.TestCase):
    def test_compose_iz_ops_combines_duplicate_terms(self) -> None:
        left = SparsePauliOp.from_list([("ZI", 1.0), ("IZ", 2.0)])
        right = SparsePauliOp.from_list([("IZ", 3.0), ("ZI", 4.0)])

        actual = compose_IZ_ops(left, right).simplify()
        expected = SparsePauliOp.from_list([("ZZ", 11.0), ("II", 10.0)]).simplify()

        self.assertTrue(actual.equiv(expected))

    def test_sampler_v2_counts_are_extracted(self) -> None:
        self.assertEqual(
            ProteinFoldingProblem._extract_counts(_PrimitiveResult()),
            {"001": 2, "111": 5},
        )

    def test_plain_count_mapping_is_accepted(self) -> None:
        counts = {"000": 1, "101": 3}
        self.assertIs(ProteinFoldingProblem._extract_counts(counts), counts)

    def test_vqec_optimizers_accept_estimator_v2(self) -> None:
        ansatz = real_amplitudes(2, reps=1)
        observable = SparsePauliOp("ZI")
        constraint = SparsePauliOp("IZ")
        estimator = StatevectorEstimator()
        gradient = ParamShiftEstimatorGradient(estimator)
        parameters = np.zeros(ansatz.num_parameters)

        for optimizer_type in (PerturbedPrimalDualOpt, OptimisticGDAOpt):
            optimizer = optimizer_type(
                observable, [constraint], ansatz, estimator, gradient
            )
            self.assertAlmostEqual(
                optimizer.get_expectation(ansatz, observable, parameters), 1.0
            )

    def test_chance_vqec_accepts_separate_sampling_circuit(self) -> None:
        class Problem:
            num_qubits = 2
            objective_scale = 1.0
            constraint_count = 1

        class Sampler:
            def run(self, *args, **kwargs):  # pragma: no cover - not submitted
                raise AssertionError("constructor test must not submit")

        ansatz = real_amplitudes(2, reps=1)
        sampling_circuit = ansatz.copy()
        sampling_circuit.measure_all()
        solver = ChanceConstrainedVQEC(
            Problem(),
            ansatz,
            Sampler(),
            constraint_limits=[0.01],
            shots=100,
            sampling_circuit=sampling_circuit,
        )
        self.assertIs(solver._measured_ansatz, sampling_circuit)


if __name__ == "__main__":
    unittest.main()
