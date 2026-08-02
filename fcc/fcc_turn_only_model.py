# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

"""Sample-wise objective and overlap indicators for compact FCC turns."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from functools import lru_cache
from typing import Protocol

import numpy as np
from qiskit.quantum_info import SparsePauliOp

from .fcc_mj_interaction import MiyazawaJerniganInteraction
from .fcc_penalty_parameters import PenaltyParameters
from .fcc_peptide import Peptide
from .fcc_protein_shape import (
    compact_turn_qubit_blocks,
    compact_turn_qubit_count,
    decode_compact_turn_bitstring,
    noncovalent_residue_pairs,
    turn_sequence_to_lattice_positions,
)
from .fcc_qubit_op_builder import QubitOpBuilder


DEFAULT_INTERACTION_MULTIPLIER = 1.0
DEFAULT_OBJECTIVE_SCALE = 1.0


class PairInteraction(Protocol):
    """Structural type for residue-pair interaction models."""

    def calculate_energy_matrix(self, residue_sequence: str) -> np.ndarray:
        """Return a residue-pair energy matrix."""


def _read_only(values: np.ndarray) -> np.ndarray:
    values.setflags(write=False)
    return values


@dataclass(frozen=True)
class TurnOnlyBitstringEvaluation:
    """Objective decomposition and indicators for one compact bitstring."""

    bitstring: str
    turn_sequence: tuple[int | None, ...]
    positions: np.ndarray
    squared_distances: np.ndarray
    contact_indicators: np.ndarray
    overlap_indicators: np.ndarray
    backtracking_energy: float
    redundancy_energy: float
    contact_energy: float
    objective: float
    physical_encoding: bool
    fully_valid: bool

    @property
    def any_overlap(self) -> bool:
        """Whether any configured residue pair overlaps."""

        return bool(np.any(self.overlap_indicators))


@dataclass(frozen=True)
class TurnOnlyBitstringBatch:
    """Vectorized evaluations for a caller-supplied bitstring collection."""

    bitstrings: tuple[str, ...]
    turn_sequences: tuple[tuple[int | None, ...], ...]
    positions: np.ndarray
    squared_distances: np.ndarray
    contact_indicators: np.ndarray
    overlap_indicators: np.ndarray
    backtracking_energies: np.ndarray
    redundancy_energies: np.ndarray
    contact_energies: np.ndarray
    objectives: np.ndarray
    physical_mask: np.ndarray
    valid_mask: np.ndarray

    @property
    def any_overlap_mask(self) -> np.ndarray:
        """Return a read-only joint-overlap indicator for each bitstring."""

        if self.overlap_indicators.shape[1] == 0:
            values = np.zeros(len(self.bitstrings), dtype=bool)
        else:
            values = np.any(self.overlap_indicators, axis=1)
        return _read_only(values)


@dataclass(eq=False, frozen=True)
class TurnOnlyFCCModel:
    """Sampled diagonal model for a fixed peptide's compact FCC register.

    Bitstrings use Qiskit's displayed most-significant-bit-first order. Unused
    four-bit turn codes receive zero displacement and remain in the sampled
    domain, where ``H_redun`` penalizes them. Samples are never postselected.
    """

    sequence: str
    num_qubits: int
    turn_qubit_blocks: tuple[tuple[int, ...], ...]
    noncovalent_pairs: tuple[tuple[int, int], ...]
    constrained_pairs: tuple[tuple[int, int], ...]
    constrained_pair_indices: tuple[int, ...]
    pair_energies: np.ndarray
    backtracking_hamiltonian: SparsePauliOp
    redundancy_hamiltonian: SparsePauliOp
    objective_scale: float = DEFAULT_OBJECTIVE_SCALE

    @property
    def peptide_length(self) -> int:
        """Number of residues in the modeled peptide."""

        return len(self.sequence)

    @property
    def constraint_count(self) -> int:
        """Number of pairwise exact-overlap chance constraints."""

        return len(self.constrained_pairs)

    def evaluate_bitstring(self, bitstring: str) -> TurnOnlyBitstringEvaluation:
        """Score one compact turn bitstring."""

        return self._evaluate_bitstring(self._validated_bitstring(bitstring))

    def objective_value(self, bitstring: str) -> float:
        """Return ``H_back + H_redun + E_int`` for one bitstring."""

        return self.evaluate_bitstring(bitstring).objective

    def contact_indicator_values(self, bitstring: str) -> np.ndarray:
        """Return ``1[D_mn == 2]`` for all objective contact pairs."""

        return self.evaluate_bitstring(bitstring).contact_indicators

    def overlap_indicator_values(self, bitstring: str) -> np.ndarray:
        """Return ``1[D_mn == 0]`` for configured chance constraints."""

        return self.evaluate_bitstring(bitstring).overlap_indicators

    def evaluate_bitstrings(
        self, bitstrings: Iterable[str]
    ) -> TurnOnlyBitstringBatch:
        """Score a finite bitstring collection, preserving order and duplicates."""

        keys = tuple(self._validated_bitstring(value) for value in bitstrings)
        if not keys:
            raise ValueError("At least one bitstring is required")
        rows = tuple(self._evaluate_bitstring(key) for key in keys)
        return TurnOnlyBitstringBatch(
            bitstrings=keys,
            turn_sequences=tuple(row.turn_sequence for row in rows),
            positions=_read_only(np.stack([row.positions for row in rows])),
            squared_distances=_read_only(
                np.stack([row.squared_distances for row in rows])
            ),
            contact_indicators=_read_only(
                np.stack([row.contact_indicators for row in rows])
            ),
            overlap_indicators=_read_only(
                np.stack([row.overlap_indicators for row in rows])
            ),
            backtracking_energies=_read_only(
                np.asarray([row.backtracking_energy for row in rows], dtype=float)
            ),
            redundancy_energies=_read_only(
                np.asarray([row.redundancy_energy for row in rows], dtype=float)
            ),
            contact_energies=_read_only(
                np.asarray([row.contact_energy for row in rows], dtype=float)
            ),
            objectives=_read_only(
                np.asarray([row.objective for row in rows], dtype=float)
            ),
            physical_mask=_read_only(
                np.asarray([row.physical_encoding for row in rows], dtype=bool)
            ),
            valid_mask=_read_only(
                np.asarray([row.fully_valid for row in rows], dtype=bool)
            ),
        )

    def bitstring_primitive_values(self, bitstrings: Iterable[str]) -> np.ndarray:
        """Return normalized objective and overlap indicators per bitstring."""

        batch = self.evaluate_bitstrings(bitstrings)
        values = np.column_stack(
            (
                batch.objectives / self.objective_scale,
                batch.overlap_indicators.astype(float),
            )
        )
        return _read_only(values)

    def primitive_values(
        self,
        counts: Mapping[str, int | float]
        | Sequence[Mapping[str, int | float]],
    ) -> np.ndarray:
        """Return sampled objective and overlap-probability expectations."""

        mappings = (counts,) if isinstance(counts, Mapping) else tuple(counts)
        if not mappings:
            raise ValueError("At least one counts mapping is required")
        expectations = []
        for sample in mappings:
            keys, weights = self._validated_counts(sample)
            per_bitstring = self.bitstring_primitive_values(keys)
            expectations.append(weights @ per_bitstring / float(weights.sum()))
        return np.asarray(expectations, dtype=float)

    def counts_distribution(
        self, counts: Mapping[str, int | float]
    ) -> tuple[tuple[str, ...], np.ndarray, float]:
        """Return canonical keys, probabilities, and total sample weight."""

        keys, weights = self._validated_counts(counts)
        total = float(weights.sum())
        return keys, _read_only(weights / total), total

    def _validated_bitstring(self, bitstring: str) -> str:
        if not isinstance(bitstring, str):
            raise TypeError("Compact turn bitstrings must be strings")
        value = bitstring.replace(" ", "")
        if len(value) != self.num_qubits:
            raise ValueError(
                f"Expected {self.num_qubits} compact turn bits, received {len(value)}"
            )
        if set(value) - {"0", "1"}:
            raise ValueError("A compact turn bitstring may contain only '0' and '1'")
        return value

    def _validated_counts(
        self, counts: Mapping[str, int | float]
    ) -> tuple[tuple[str, ...], np.ndarray]:
        if not isinstance(counts, Mapping) or not counts:
            raise ValueError("Counts must be a nonempty mapping")
        merged: dict[str, float] = {}
        for key, count in counts.items():
            bitstring = self._validated_bitstring(key)
            value = float(count)
            if not math.isfinite(value) or value < 0.0:
                raise ValueError("Counts must be finite and nonnegative")
            merged[bitstring] = merged.get(bitstring, 0.0) + value
        merged = {key: value for key, value in merged.items() if value > 0.0}
        if not merged:
            raise ValueError("Counts must contain positive total weight")
        keys = tuple(merged)
        return keys, np.asarray([merged[key] for key in keys], dtype=float)

    @lru_cache(maxsize=32768)
    def _evaluate_bitstring(self, bitstring: str) -> TurnOnlyBitstringEvaluation:
        turns = decode_compact_turn_bitstring(bitstring, self.peptide_length)
        positions = turn_sequence_to_lattice_positions(
            turns, invalid_turns_as_zero=True
        )
        distances = np.empty(len(self.noncovalent_pairs), dtype=np.int64)
        for index, (lower, upper) in enumerate(self.noncovalent_pairs):
            difference = positions[lower] - positions[upper]
            distances[index] = int(difference @ difference)

        contacts = distances == 2
        contact_energy = float(self.pair_energies @ contacts.astype(float))
        overlap = (
            distances[np.asarray(self.constrained_pair_indices, dtype=int)] == 0
            if self.constrained_pair_indices
            else np.empty(0, dtype=bool)
        )
        backtracking = _diagonal_pauli_value(
            self.backtracking_hamiltonian, bitstring
        )
        redundancy = _diagonal_pauli_value(self.redundancy_hamiltonian, bitstring)
        physical = all(turn is not None for turn in turns)
        self_avoiding = len({tuple(position) for position in positions}) == len(
            positions
        )
        objective = backtracking + redundancy + contact_energy

        return TurnOnlyBitstringEvaluation(
            bitstring=bitstring,
            turn_sequence=turns,
            positions=_read_only(np.asarray(positions, dtype=np.int16)),
            squared_distances=_read_only(distances),
            contact_indicators=_read_only(contacts),
            overlap_indicators=_read_only(np.asarray(overlap, dtype=bool)),
            backtracking_energy=float(backtracking),
            redundancy_energy=float(redundancy),
            contact_energy=contact_energy,
            objective=float(objective),
            physical_encoding=physical,
            fully_valid=physical and self_avoiding,
        )


def _diagonal_pauli_value(operator: SparsePauliOp, bitstring: str) -> float:
    """Evaluate an I/Z operator on a basis bitstring without a statevector."""

    if operator.num_qubits != len(bitstring):
        raise ValueError("Hamiltonian and bitstring qubit counts differ")
    if np.any(operator.paulis.x):
        raise ValueError("Turn-penalty Hamiltonians must be diagonal in I/Z")
    little_endian_bits = np.fromiter(
        (character == "1" for character in reversed(bitstring)), dtype=np.uint8
    )
    parities = (
        operator.paulis.z.astype(np.uint8) @ little_endian_bits
    ) % np.uint8(2)
    eigenvalues = 1.0 - 2.0 * parities.astype(float)
    value = np.dot(operator.coeffs, eigenvalues)
    if abs(float(np.imag(value))) > 1e-10:
        raise ValueError("Diagonal Hamiltonian produced a complex expectation")
    return float(np.real(value))


def _sequence_value(peptide: Peptide | str) -> str:
    sequence = (
        peptide.peptide_sequence if isinstance(peptide, Peptide) else str(peptide)
    )
    if len(sequence) < 3:
        raise ValueError("The compact FCC encoding requires at least three residues")
    return sequence


def _validated_pairs(
    pairs: Iterable[tuple[int, int]], peptide_length: int
) -> tuple[tuple[int, int], ...]:
    result = tuple((int(lower), int(upper)) for lower, upper in pairs)
    if len(set(result)) != len(result):
        raise ValueError("Residue-pair list contains duplicates")
    if any(
        lower < 0
        or upper >= peptide_length
        or lower >= upper
        or upper - lower < 2
        for lower, upper in result
    ):
        raise ValueError("Residue pairs must be noncovalent and in ascending order")
    return result


def build_turn_only_fcc_model(
    peptide: Peptide | str,
    *,
    interaction: PairInteraction | None = None,
    penalty_parameters: PenaltyParameters | None = None,
    interaction_multiplier: float = DEFAULT_INTERACTION_MULTIPLIER,
    constrained_pairs: Iterable[tuple[int, int]] | None = None,
    objective_scale: float = DEFAULT_OBJECTIVE_SCALE,
) -> TurnOnlyFCCModel:
    """Build a turn-only sample scorer for a peptide.

    The raw objective is ``H_back + H_redun + sum(epsilon_mn * 1[D_mn == 2])``.
    Chance primitives are the distinct indicators ``1[D_mn == 0]``. No
    contact ancilla, contact Hamiltonian, statevector, or basis enumeration is
    constructed.
    """

    sequence = _sequence_value(peptide)
    peptide_object = peptide if isinstance(peptide, Peptide) else Peptide(sequence)
    peptide_length = len(sequence)
    multiplier = float(interaction_multiplier)
    scale = float(objective_scale)
    if not math.isfinite(multiplier):
        raise ValueError("interaction_multiplier must be finite")
    if not math.isfinite(scale) or scale <= 0.0:
        raise ValueError("objective_scale must be finite and positive")

    penalties = penalty_parameters or PenaltyParameters()
    if not all(
        math.isfinite(float(value)) and float(value) >= 0.0
        for value in (penalties.penalty_back, penalties.penalty_redun)
    ):
        raise ValueError("Turn penalties must be finite and nonnegative")

    noncovalent_pairs = noncovalent_residue_pairs(
        peptide_length, minimum_separation=2
    )
    default_constraints = noncovalent_residue_pairs(
        peptide_length, minimum_separation=3
    )
    constraints = _validated_pairs(
        default_constraints if constrained_pairs is None else constrained_pairs,
        peptide_length,
    )
    if not set(constraints).issubset(noncovalent_pairs):
        raise ValueError("Constrained pairs must be a subset of contact pairs")
    pair_indices = {pair: index for index, pair in enumerate(noncovalent_pairs)}
    constrained_pair_indices = tuple(pair_indices[pair] for pair in constraints)

    interaction_model = interaction or MiyazawaJerniganInteraction("mj_matrix")
    energy_matrix = np.asarray(
        interaction_model.calculate_energy_matrix(sequence), dtype=float
    )
    if energy_matrix.shape != (peptide_length, peptide_length):
        raise ValueError("Interaction matrix shape does not match the peptide")
    pair_energies = multiplier * np.asarray(
        [energy_matrix[pair] for pair in noncovalent_pairs], dtype=float
    )
    if not np.isfinite(pair_energies).all():
        raise ValueError("Pair interaction coefficients must be finite")

    builder = QubitOpBuilder(
        peptide_object,
        energy_matrix,
        penalties,
        build_geometry_maps=False,
    )
    h_back, h_redun = builder.build_turn_penalty_ops()
    num_qubits = compact_turn_qubit_count(peptide_length)
    if h_back.num_qubits != num_qubits or h_redun.num_qubits != num_qubits:
        raise RuntimeError("Compact penalty Hamiltonian has the wrong qubit count")

    return TurnOnlyFCCModel(
        sequence=sequence,
        num_qubits=num_qubits,
        turn_qubit_blocks=compact_turn_qubit_blocks(peptide_length),
        noncovalent_pairs=noncovalent_pairs,
        constrained_pairs=constraints,
        constrained_pair_indices=constrained_pair_indices,
        pair_energies=_read_only(pair_energies),
        backtracking_hamiltonian=h_back,
        redundancy_hamiltonian=h_redun,
        objective_scale=scale,
    )
