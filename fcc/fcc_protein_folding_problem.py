# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

"""A class for the protein folding problem on an FCC lattice."""

from qiskit.quantum_info import SparsePauliOp
from qiskit_algorithms.minimum_eigensolvers import SamplingVQEResult
from .fcc_peptide import Peptide
from .fcc_penalty_parameters import PenaltyParameters
from .fcc_qubit_op_builder import QubitOpBuilder
from .fcc_protein_folding_result import ProteinFoldingResult
from .measurement_utils import (
    extract_counts_from_result,
    normalize_bitstring_mapping,
    normalize_counts,
)
from .utils import remove_unused_qubits
from typing import Any, Protocol


class InteractionModel(Protocol):
    """Protocol for interaction models that can generate pair-energy matrices."""

    def calculate_energy_matrix(self, residue_sequence: str): ...

    def validate_residue_sequence(self, residue_sequence: str) -> None: ...


class ProteinFoldingProblem:

    def __init__(
        self,
        peptide: Peptide,
        interaction: InteractionModel,
        penalty_parameters: PenaltyParameters,
    ):
        """
        Args:
            peptide: A Peptide object that includes all information about a
            protein.
            interaction: An interaction model that defines the energy matrix
            through ``calculate_energy_matrix``.
            penalty_parameters: Parameters that define the strength of
            constraints in the problem.
        """
        self._peptide = peptide
        self._interaction = interaction
        validator = getattr(interaction, "validate_residue_sequence", None)
        if callable(validator):
            validator(peptide.peptide_sequence)
        self._penalty_parameters = penalty_parameters
        self._pair_energies = interaction.calculate_energy_matrix(
            peptide.peptide_sequence
        )
        self._qubit_op_builder = QubitOpBuilder(
            self._peptide,
            self._pair_energies,
            self._penalty_parameters,
        )
        self._unused_qubits: list[int] | None = None
        self._reduced_qubit_op: SparsePauliOp | None = None

    def qubit_op(self, r2_threshold: float = 1.0, chunk: int = 20) -> SparsePauliOp:
        """
        Builds the total qubit operator for the full Hamiltonian encoding a
        protein folding problem on the FCC lattice.

        Args:
            r2_threshold: The threshold for the R^2 score of the Chebyshev fit
            when building the non-overlapping constraint.
            chunk: Size of the chunks to split the qubit operator into.

        Returns:
            A qubit operator for the full Hamiltonian encoding a protein folding
            problem.
        """
        qubit_op = self._qubit_op_builder.build_qubit_op(r2_threshold, chunk)
        reduced_qubit_op, unused_qubits = remove_unused_qubits(qubit_op)
        self._unused_qubits = unused_qubits
        self._reduced_qubit_op = reduced_qubit_op
        return reduced_qubit_op

    def olap_constr_ops(self) -> tuple[list[tuple[int, int]], list[SparsePauliOp]]:
        """
        Builds the overlap constraint operators for the protein folding problem
        on the FCC lattice that are used in the VQEC approach.

        This method is only available when ``penalty_olap`` is ``None``. When a
        finite overlap penalty is baked into the Hamiltonian, there is no
        separate set of dualized overlap constraints to return.

        Call :meth:`qubit_op` first so the overlap constraints are compressed
        using the exact same unused-qubit map as the Hamiltonian passed to the
        solver.

        Returns:
            A tuple of a list of bead pairs and a list of corresponding qubit
            operators for the overlap constraints.

        Raises:
            ValueError: If ``penalty_olap`` is not ``None``.
        """
        if self._penalty_parameters.penalty_olap is not None:
            raise ValueError(
                "olap_constr_ops() is only available when penalty_olap is None. "
                "Set penalty_olap=None to build explicit overlap constraints for VQEC."
            )
        if self._unused_qubits is None:
            raise RuntimeError(
                "Call qubit_op(...) before olap_constr_ops() so the constraint operators use the same qubit compression as the Hamiltonian."
            )

        olap_constr_ops_dict = self._qubit_op_builder.build_olap_constr_ops()
        for pair in olap_constr_ops_dict:
            olap_constr_ops_dict[pair], _ = remove_unused_qubits(
                olap_constr_ops_dict[pair], self._unused_qubits
            )
        olap_constr_ops = list(olap_constr_ops_dict.values())
        bead_pairs = list(olap_constr_ops_dict.keys())
        return bead_pairs, olap_constr_ops

    def interpret(self, raw_result: SamplingVQEResult | Any) -> ProteinFoldingResult:
        """
        Interprets the raw algorithm result and returns a ProteinFoldingResult
        object.

        Args:
            raw_result: A raw result of the protein folding problem.

        Returns:
            A ProteinFoldingResult object that includes the interpreted result.
        """
        if self._unused_qubits is None or self._reduced_qubit_op is None:
            raise RuntimeError(
                "Call qubit_op(...) before interpret() so result decoding uses the same compressed Hamiltonian that produced the raw result."
            )
        num_active_qubits = self._reduced_qubit_op.num_qubits
        try:
            best_turn_bitstring = next(
                iter(
                    normalize_counts(
                        {raw_result.best_measurement["bitstring"]: 1}, num_active_qubits
                    )
                )
            )
        except (AttributeError, KeyError, TypeError):
            if hasattr(raw_result, "quasi_dists"):
                prob_dist = normalize_bitstring_mapping(
                    raw_result.quasi_dists[0].binary_probabilities(), num_active_qubits
                )
                # Find the most probable bitstring
                best_turn_bitstring = max(prob_dist, key=prob_dist.get)
            else:
                try:
                    counts_like = extract_counts_from_result(raw_result)
                except Exception:
                    try:
                        counts_like = extract_counts_from_result(raw_result[0])
                    except Exception as exc:
                        raise TypeError(
                            "Unsupported result type. Expected a SamplingVQEResult, a legacy SamplerResult, or a SamplerV2 PrimitiveResult."
                        ) from exc
                counts = normalize_counts(counts_like, num_active_qubits)
                best_turn_bitstring = max(counts, key=counts.get)
        return ProteinFoldingResult(
            peptide=self._peptide,
            unused_qubits=self.unused_qubits,
            solution_bitstring=best_turn_bitstring,
        )

    @property
    def unused_qubits(self) -> list[int]:
        """Returns the list of indices for qubits in the original problem
        formulation that were removed during compression."""
        return [] if self._unused_qubits is None else self._unused_qubits

    @property
    def peptide(self) -> Peptide:
        """Returns the peptide defining the protein subject to the folding
        problem."""
        return self._peptide
