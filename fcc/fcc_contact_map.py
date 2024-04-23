"""A class for building a contact map that stores contacts between beads in a peptide on the FCC lattice."""

from collections import defaultdict
from qiskit.quantum_info import SparsePauliOp
from fcc_peptide import Peptide
from utils import build_full_identity, build_pauli_z_op


class ContactMap:

    def __init__(self, peptide: Peptide):
        """
        Args:
            peptide: A Peptide object that includes all information about a protein.
        """
        self._peptide = peptide
        self._peptide_length = peptide.peptide_length
        self._num_qubits = pow(self._peptide_length - 1, 2)
        (self._contact_map, self._num_contacts) = self._build_contact_map(peptide)

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

    def _build_contact_map(
        self, peptide: Peptide
    ) -> tuple[defaultdict[int, dict[int, SparsePauliOp]], int]:
        """
        Builds a contact map for a given peptide -- a list of Pauli operators that represent nearest neighbor interactions. A nearest neighbor interaction between 2 beads is encoded using 1 qubit. Since we are considering only the main chain, the maximum number of qubits needed to encode all interactions is :math:`(N-1)^2`, where :math:`N` is the number of beads in the main chain. We build the contact operators based on the following indexing: lower_bead_idx * (N-1) + upper_bead_idx.
        Note that some qubits may not be necessary because of no contacts between respective beads. Such qubits are deleted in the compression phase happening during the creation of a final qubit operator for the problem.

        Args:
            peptide: A Peptide object that includes all information about a protein.

        Returns:
            A tuple of a contact map and the number of contacts calculated.
        """
        num_contacts = 0
        contact_map = defaultdict(dict)
        for lower_bead_idx in range(self._peptide_length - 2):
            for upper_bead_idx in range(lower_bead_idx + 2, self._peptide_length):
                contact_map[lower_bead_idx][upper_bead_idx] = self._create_contact_op(
                    lower_bead_idx, upper_bead_idx, self._num_qubits
                )
                num_contacts += 1
        return contact_map, num_contacts

    def _create_contact_op(
        self, lower_bead_idx: int, upper_bead_idx: int, num_qubits: int
    ) -> SparsePauliOp:
        """
        Creates a Pauli operator for a nearest neighbor interaction between two beads.

        Args:
            lower_bead_idx: The index of the lower bead.
            upper_bead_idx: The index of the upper bead.
            num_qubits: The number of qubits needed to encode all interactions.

        Returns:
            A Pauli operator for a nearest neighbor interaction between two beads.
        """
        z_op_idx = lower_bead_idx * (self._peptide_length - 1) + upper_bead_idx
        contact_op = (
            (build_full_identity(num_qubits) - build_pauli_z_op(num_qubits, {z_op_idx}))
            / 2
        ).simplify()
        return contact_op
