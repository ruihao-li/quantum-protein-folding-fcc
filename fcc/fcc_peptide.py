# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

"""A class defining the main chain of a peptide."""

from qiskit.quantum_info import SparsePauliOp
from .fcc_bead import Bead
from .utils import build_full_identity, build_pauli_z_op


class Peptide:

    def __init__(self, peptide_sequence: str):
        """
        Args:
            peptide_sequence: String of characters that define residues for the
            main chain of a peptide. Valid residue types are [A, C, D, E, F, G,
            H, I, K, L, M, N, P, Q, R, S, T, V, W, Y].
        """
        self._peptide_sequence = peptide_sequence
        self._beads_list = self._build_main_chain(peptide_sequence)

    @property
    def peptide_sequence(self) -> str:
        """Returns a residue sequence for the peptide chain."""
        return self._peptide_sequence

    @property
    def beads_list(self) -> list[Bead]:
        """Returns the list of all beads in the chain."""
        return self._beads_list

    @property
    def peptide_length(self) -> int:
        """Returns the length of the peptide."""
        return len(self._peptide_sequence)

    def _build_main_chain(self, peptide_sequence: str) -> list[Bead]:
        """
        Creates a main chain for a given peptide sequence.

        Args:
            peptide_sequence: String of characters that define residues for the
            main chain of a peptide.

        Returns:
            A list of MainBead instances.
        """
        main_chain = []
        main_chain_len = len(peptide_sequence)
        for main_bead_id in range(main_chain_len - 1):
            bead_turn_qubit_0 = self._build_turn_qubit(main_chain_len, 4 * main_bead_id)
            bead_turn_qubit_1 = self._build_turn_qubit(
                main_chain_len, 4 * main_bead_id + 1
            )
            bead_turn_qubit_2 = self._build_turn_qubit(
                main_chain_len, 4 * main_bead_id + 2
            )
            bead_turn_qubit_3 = self._build_turn_qubit(
                main_chain_len, 4 * main_bead_id + 3
            )
            main_bead = Bead(
                main_bead_id,
                peptide_sequence[main_bead_id],
                (
                    bead_turn_qubit_0,
                    bead_turn_qubit_1,
                    bead_turn_qubit_2,
                    bead_turn_qubit_3,
                ),
            )
            main_chain.append(main_bead)
        # Add the last bead
        last_bead = Bead(main_chain_len - 1, peptide_sequence[-1], None)
        main_chain.append(last_bead)
        return main_chain

    @staticmethod
    def _build_turn_qubit(chain_len: int, pauli_z_index: int) -> SparsePauliOp:
        """
        Builds a SparsePauliOp of length 4 * (chain_len - 1) (number of qubits
        necessary to encode all turns for the chain of length chain_len on an
        FCC lattice) with a Pauli Z operator at a given index: :math:`q_i = (I -
        Z_i)/2`.

        Args:
            chain_len: Length of the chain.
            pauli_z_index: Index of a Pauli Z operator in a turn operator.

        Returns:
            A Pauli operator that encodes the turn following from a given bead
            index.
        """
        num_turn_qubits = 4 * (chain_len - 1)
        norm_factor = 0.5
        turn_qubit = (
            norm_factor * build_full_identity(num_turn_qubits)
            - norm_factor * build_pauli_z_op(num_turn_qubits, {pauli_z_index})
        ).simplify()
        return turn_qubit
