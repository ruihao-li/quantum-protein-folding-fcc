"""Characterize the historical full-register/PolyFit path."""

from __future__ import annotations

import unittest

import numpy as np

from fcc import (
    MiyazawaJerniganInteraction,
    Peptide,
    PenaltyParameters,
    ProteinFoldingProblem,
)


class PolyFitCompatibilityTests(unittest.TestCase):
    def test_default_builder_remains_eager_and_full_register(self) -> None:
        problem = ProteinFoldingProblem(
            peptide=Peptide("AAA"),
            interaction=MiyazawaJerniganInteraction("mj_matrix"),
            penalty_parameters=PenaltyParameters(
                penalty_back=10.0,
                penalty_redun=10.0,
                penalty_olap=None,
            ),
        )
        builder = problem._qubit_op_builder
        self.assertIsNotNone(builder._contact_map)
        self.assertIsNotNone(builder._distance_map)
        self.assertEqual(builder._contact_map.num_contacts, 1)
        self.assertEqual(builder._distance_map.num_distances, 3)

        operator = problem.qubit_op()
        self.assertEqual(operator.num_qubits, 3)
        self.assertEqual(operator.size, 6)
        self.assertEqual(problem.unused_qubits, [0, 1, 2, 3, 6, 7])
        coefficients = {label: value.real for label, value in operator.to_list()}
        expected = {
            "III": 2.72,
            "IZI": 2.72,
            "IZZ": 1.36,
            "ZII": -2.72,
            "ZZI": -2.72,
            "ZZZ": -1.36,
        }
        self.assertEqual(set(coefficients), set(expected))
        for label, value in expected.items():
            self.assertTrue(np.isclose(coefficients[label], value))

    def test_legacy_pairwise_distance_constraint_is_unchanged(self) -> None:
        problem = ProteinFoldingProblem(
            peptide=Peptide("ACDE"),
            interaction=MiyazawaJerniganInteraction("mj_matrix"),
            penalty_parameters=PenaltyParameters(
                penalty_back=10.0,
                penalty_redun=10.0,
                penalty_olap=None,
            ),
        )
        operator = problem.qubit_op()
        pairs, constraints = problem.olap_constr_ops()
        self.assertEqual(operator.num_qubits, 9)
        self.assertEqual(operator.size, 107)
        self.assertEqual(pairs, [(0, 3)])
        self.assertEqual(len(constraints), 1)
        self.assertEqual(constraints[0].num_qubits, 9)
        self.assertEqual(constraints[0].size, 35)


if __name__ == "__main__":
    unittest.main()
