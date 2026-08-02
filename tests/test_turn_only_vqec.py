"""Deterministic regression tests for compact turn-only FCC scoring."""

from __future__ import annotations

import unittest
from unittest.mock import patch

import numpy as np
from qiskit.quantum_info import Statevector

from fcc import (
    Peptide,
    PenaltyParameters,
    ProteinShapeDecoder,
    ProteinShapeFileGen,
    build_turn_only_fcc_model,
    compact_turn_qubit_count,
    decode_compact_turn_bitstring,
    noncovalent_residue_pairs,
    turn_sequence_to_lattice_positions,
)
from fcc.fcc_protein_shape import turn_sequence_to_compact_bitstring


class ConstantInteraction:
    """Return a caller-selected coefficient for every residue pair."""

    def __init__(self, coefficient: float = 1.0) -> None:
        self.coefficient = coefficient

    def calculate_energy_matrix(self, residue_sequence: str) -> np.ndarray:
        size = len(residue_sequence)
        return np.full((size, size), self.coefficient, dtype=float)


class TurnOnlyFCCModelTests(unittest.TestCase):
    def test_compact_register_and_general_constraint_pairs(self) -> None:
        for length in (3, 4, 5, 6, 8):
            self.assertEqual(compact_turn_qubit_count(length), 4 * length - 10)
        model = build_turn_only_fcc_model("KLVFFA")
        self.assertEqual(model.num_qubits, 14)
        self.assertEqual(
            model.constrained_pairs,
            ((0, 3), (0, 4), (0, 5), (1, 4), (1, 5), (2, 5)),
        )
        self.assertEqual(
            model.constrained_pairs,
            noncovalent_residue_pairs(6, minimum_separation=3),
        )

    def test_turn_only_build_never_constructs_geometry_maps(self) -> None:
        with patch(
            "fcc.fcc_qubit_op_builder.ContactMap",
            side_effect=AssertionError("contact map constructed"),
        ), patch(
            "fcc.fcc_qubit_op_builder.DistanceMap",
            side_effect=AssertionError("distance map constructed"),
        ):
            model = build_turn_only_fcc_model("ACDE")
        self.assertEqual(model.num_qubits, 6)
        self.assertEqual(model.backtracking_hamiltonian.num_qubits, 6)
        self.assertEqual(model.redundancy_hamiltonian.num_qubits, 6)

    def test_contact_and_overlap_indicators_are_distinct(self) -> None:
        model = build_turn_only_fcc_model(
            "ACDEA",
            interaction=ConstantInteraction(),
            penalty_parameters=PenaltyParameters(3.0, 5.0),
        )
        evaluation = model.evaluate_bitstring("0000011010")

        np.testing.assert_array_equal(
            evaluation.contact_indicators, evaluation.squared_distances == 2
        )
        constrained_indices = [
            model.noncovalent_pairs.index(pair) for pair in model.constrained_pairs
        ]
        np.testing.assert_array_equal(
            evaluation.overlap_indicators,
            evaluation.squared_distances[constrained_indices] == 0,
        )
        self.assertEqual(int(evaluation.contact_indicators.sum()), 4)
        self.assertEqual(int(evaluation.overlap_indicators.sum()), 2)
        self.assertEqual(evaluation.contact_energy, 4.0)
        self.assertEqual(evaluation.objective, 4.0)

        for pair_index, distance in enumerate(evaluation.squared_distances):
            if distance == 0:
                self.assertFalse(evaluation.contact_indicators[pair_index])

    def test_objective_includes_backtracking_and_redundancy(self) -> None:
        penalties = PenaltyParameters(3.0, 5.0)
        interaction = ConstantInteraction()

        redundant_model = build_turn_only_fcc_model(
            "ACDEA", interaction=interaction, penalty_parameters=penalties
        )
        redundant = redundant_model.evaluate_bitstring("0001010010")
        self.assertFalse(redundant.physical_encoding)
        self.assertEqual(redundant.redundancy_energy, 5.0)
        self.assertAlmostEqual(
            redundant.objective,
            redundant.backtracking_energy
            + redundant.redundancy_energy
            + redundant.contact_energy,
        )

        backtracking_model = build_turn_only_fcc_model(
            "KLVFFA", interaction=interaction, penalty_parameters=penalties
        )
        backtracking = backtracking_model.evaluate_bitstring("11111111111111")
        self.assertEqual(backtracking.backtracking_energy, 3.0)
        self.assertAlmostEqual(
            backtracking.objective,
            backtracking.backtracking_energy
            + backtracking.redundancy_energy
            + backtracking.contact_energy,
        )

    def test_penalty_values_match_basis_state_expectations(self) -> None:
        model = build_turn_only_fcc_model("KLVFFA")
        for bitstring in ("00000000000000", "00000000010001", "11111111111111"):
            evaluation = model.evaluate_bitstring(bitstring)
            state = Statevector.from_label(bitstring)
            self.assertAlmostEqual(
                evaluation.backtracking_energy,
                float(
                    np.real(
                        state.expectation_value(model.backtracking_hamiltonian)
                    )
                ),
            )
            self.assertAlmostEqual(
                evaluation.redundancy_energy,
                float(
                    np.real(state.expectation_value(model.redundancy_hamiltonian))
                ),
            )

    def test_counts_aggregation_keeps_invalid_samples(self) -> None:
        model = build_turn_only_fcc_model("ACDEA", objective_scale=10.0)
        counts = {"0000011010": 3, "0001010010": 1}
        per_bitstring = model.bitstring_primitive_values(counts)
        expected = np.asarray([0.75, 0.25]) @ per_bitstring
        np.testing.assert_allclose(
            model.primitive_values(counts)[0], expected, atol=1e-14, rtol=0.0
        )
        self.assertFalse(model.evaluate_bitstring("0001010010").physical_encoding)

    def test_legacy_decoder_and_coordinates_share_compact_geometry(self) -> None:
        peptide = Peptide("ACDE")
        turns = (0, 8, 4)
        bitstring = turn_sequence_to_compact_bitstring(turns)
        self.assertEqual(decode_compact_turn_bitstring(bitstring, 4), turns)
        decoder = ProteinShapeDecoder(peptide, bitstring)
        self.assertEqual(tuple(decoder.turn_sequence), turns)
        generated = ProteinShapeFileGen(
            peptide, list(turns)
        ).amino_acid_positions
        expected = 3.8 / np.sqrt(2.0) * turn_sequence_to_lattice_positions(turns)
        np.testing.assert_allclose(generated, expected, atol=0.0, rtol=0.0)


if __name__ == "__main__":
    unittest.main()
