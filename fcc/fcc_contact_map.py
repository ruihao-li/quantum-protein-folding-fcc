# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

"""A class for building a contact map that stores contacts between beads in a peptide on the FCC lattice."""

from collections import defaultdict
from qiskit.quantum_info import SparsePauliOp
from .fcc_peptide import Peptide
from .utils import build_full_identity, build_pauli_z_op


class ContactMap:

    def __init__(self, peptide: Peptide):
        """
        Args:
            peptide: A Peptide object that includes all information about a protein.
        """
        self._peptide = peptide
        self._peptide_length = peptide.peptide_length
        # Only pairs of beads that are at least 2 positions apart can have contacts on the FCC lattice
        self._num_qubits = (self._peptide_length**2 - 3 * self._peptide_length + 2) // 2
        (self._contact_map, self._num_contacts) = self._build_contact_map()

    @property
    def peptide(self) -> Peptide:
        """Returns a peptide."""
        return self._peptide

    @property
    def contact_map(self) -> defaultdict[int, dict[int, SparsePauliOp]]:
        """Returns a contact map."""
        return self._contact_map

    @property
    def num_contacts(self) -> int:
        """Returns the number of contacts calculated."""
        return self._num_contacts

    @property
    def num_qubits(self) -> int:
        """Returns the number of interaction qubits in the contact map."""
        return self._num_qubits

    def _build_contact_map(
        self,
    ) -> tuple[defaultdict[int, dict[int, SparsePauliOp]], int]:
        """
        Builds a contact map for a given peptide -- a list of Pauli operators
        that represent nearest neighbor interactions. A nearest neighbor
        interaction between 2 beads is encoded using 1 qubit. Since only pairs
        of beads that are at least 2 positions apart can have contacts on the
        FCC lattice, the number of possible contacts on the main chain is
        :math:`(N^2 - 3N + 2) / 2`, where :math:`N` is the number of beads in
        the main chain. We build the contact operators based on the following
        indexing: the first pair corresponds to the last interaction qubit, the
        second pair corresponds to the second-to-last interaction qubit, and so
        on.

        Args:
            peptide: A Peptide object that includes all information about a
            protein.

        Returns:
            A tuple of a contact map and the number of contacts calculated.
        """
        num_contacts = 0
        contact_map = defaultdict(dict)
        for lower_bead_idx in range(self._peptide_length - 2):
            for upper_bead_idx in range(lower_bead_idx + 2, self._peptide_length):
                # pauli_z_ops: I...IZ, I...IZI, ..., ZI...I; the last qubit is the first qubit in the interaction qubit string
                contact_op = (
                    (
                        build_full_identity(self._num_qubits)
                        - build_pauli_z_op(self._num_qubits, {num_contacts})
                    )
                    / 2
                ).simplify()
                contact_map[lower_bead_idx][upper_bead_idx] = contact_op
                num_contacts += 1
        return contact_map, num_contacts
