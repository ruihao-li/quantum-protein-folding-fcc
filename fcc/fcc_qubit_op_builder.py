"""Builds the total qubit operator for the full Hamiltonian encoding a protein folding problem on the FCC lattice."""

import numpy as np
from qiskit.quantum_info import SparsePauliOp
from fcc_peptide import Peptide
from fcc_distance_map import DistanceMap
from fcc_contact_map import ContactMap
from fcc_penalty_parameters import PenaltyParameters
from utils import (
    build_full_identity,
    build_pauli_z_op,
    fix_qubits,
)


class QubitOpBuilder:

    def __init__(
        self,
        peptide: Peptide,
        pair_energies: np.ndarray,
        penalty_parameters: PenaltyParameters,
    ):
        """
        Args:
            peptide: A Peptide object that includes all information about a protein.
            pair_energies: A matrix of pair energies between beads.
            penalty_parameters: Parameters that define the strength of constraints enforcing in the problem.
        """
        self._peptide = peptide
        self._peptide_length = peptide.peptide_length
        self._pair_energies = pair_energies
        self._penalty_parameters = penalty_parameters
        self._contact_map = ContactMap(peptide)
        self._distance_map = DistanceMap(peptide)
        self._num_config_qubits = self._distance_map._num_qubits
        self._num_contact_qubits = self._contact_map._num_qubits

    def build_qubit_op(self) -> SparsePauliOp:
        """
        Builds the total qubit operator for the full Hamiltonian encoding a protein folding problem. H_total = H_back + H_redun + H_olap + H_contact.

        Returns:
            A qubit operator for the full Hamiltonian encoding a protein folding problem.
        """
        contact_id = build_full_identity(self._num_contact_qubits)
        h_back = self._create_h_back()
        if h_back != 0:
            h_back = contact_id ^ h_back
        h_redun = self._create_h_redun()
        if h_redun != 0:
            h_redun = contact_id ^ h_redun
        h_olap = self._create_h_olap()
        if h_olap != 0:
            h_olap = contact_id ^ h_olap
        h_contact = self._create_h_contact()
        h_total = h_back + h_redun + h_olap + h_contact
        return h_total.simplify()

    def _create_turn_operator(
        self, lower_bead_idx: int, upper_bead_idx: int
    ) -> SparsePauliOp:
        """
        Creates a qubit operator for consecutive turns.

        Args:
            lower_bead_idx: The index of the lower bead.
            upper_bead_idx: The index of the upper bead.
        """
        lower_bead = self._peptide.beads_list[lower_bead_idx]
        upper_bead = self._peptide.beads_list[upper_bead_idx]

        (
            lower_turn_indicator_pxpy,
            lower_turn_indicator_pxmy,
            lower_turn_indicator_mxpy,
            lower_turn_indicator_mxmy,
            lower_turn_indicator_pypz,
            lower_turn_indicator_pymz,
            lower_turn_indicator_mypz,
            lower_turn_indicator_mymz,
            lower_turn_indicator_pxpz,
            lower_turn_indicator_pxmz,
            lower_turn_indicator_mxpz,
            lower_turn_indicator_mxmz,
        ) = lower_bead.physical_turn_indicators

        (
            upper_turn_indicator_pxpy,
            upper_turn_indicator_pxmy,
            upper_turn_indicator_mxpy,
            upper_turn_indicator_mxmy,
            upper_turn_indicator_pypz,
            upper_turn_indicator_pymz,
            upper_turn_indicator_mypz,
            upper_turn_indicator_mymz,
            upper_turn_indicator_pxpz,
            upper_turn_indicator_pxmz,
            upper_turn_indicator_mxpz,
            upper_turn_indicator_mxmz,
        ) = upper_bead.physical_turn_indicators

        turn_operator = (
            lower_turn_indicator_pxpy @ upper_turn_indicator_mxmy
            + lower_turn_indicator_pxmy @ upper_turn_indicator_mxpy
            + lower_turn_indicator_mxpy @ upper_turn_indicator_pxmy
            + lower_turn_indicator_mxmy @ upper_turn_indicator_pxpy
            + lower_turn_indicator_pypz @ upper_turn_indicator_mymz
            + lower_turn_indicator_pymz @ upper_turn_indicator_mypz
            + lower_turn_indicator_mypz @ upper_turn_indicator_pymz
            + lower_turn_indicator_mymz @ upper_turn_indicator_pypz
            + lower_turn_indicator_pxpz @ upper_turn_indicator_mxmz
            + lower_turn_indicator_pxmz @ upper_turn_indicator_mxpz
            + lower_turn_indicator_mxpz @ upper_turn_indicator_pxmz
            + lower_turn_indicator_mxmz @ upper_turn_indicator_pxpz
        ).simplify()
        return fix_qubits(turn_operator)

    def _create_h_back(self) -> SparsePauliOp:
        """
        Creates qubit operators for the H_back term, which penalizes consecutive turns in opposite directions. Note that the first bead (index 0) is omitted because the two turns following it cannot be in opposite directions by construction.

        Returns:
            A qubit operator for the H_back term.
        """
        penalty_back = self._penalty_parameters.penalty_back
        h_back = 0
        for i in range(1, self._peptide_length - 2):
            h_back += penalty_back * self._create_turn_operator(i, i + 1)
        return fix_qubits(h_back)

    def _create_h_redun(self) -> SparsePauliOp:
        """
        Creates qubit operators for the H_redun term, which penalizes redundant turns that do not correspond to any physical turns. Note that the first two turns are omitted because they are always physical by construction.

        Returns:
            A qubit operator for the H_redun term.
        """
        penalty_redun = self._penalty_parameters.penalty_redun
        h_redun = 0
        for i in range(2, self._peptide_length - 1):
            (
                turn_indicator_0010,
                turn_indicator_0001,
                turn_indicator_1101,
                turn_indicator_1110,
            ) = self._peptide.beads_list[i].unused_turn_indicators
            h_redun += penalty_redun * (
                turn_indicator_0010
                + turn_indicator_0001
                + turn_indicator_1101
                + turn_indicator_1110
            )
        return fix_qubits(h_redun)

    def _create_h_olap(self) -> SparsePauliOp:
        r"""
        Creates qubit operators for the H_olap term, which penalizes overlapping beads. To ensure the non-overlapping condition, we impose constraints on the distance function that encodes the squared distance between beads, that is, for any two beads :math:`i` and :math:`j`, we require :math:`2 \leq D_{ij} leq 2(i-j)^2`.
        Here we use the unbalanced penalization method proposed in arXiv:2211.13914 to implement the inequality constraints. The RHS constraint translates to a penalty term :math:`p_{ij}^{(1)} = -\lambda_1[2(i-j)^2-D_{ij}] + \lambda_2[2(i-j)^2-D_{ij}]^2`, where :math:`\lambda_1` and :math:`\lambda_2` are positive penalty parameters. The LHS constraint translates to a penalty term :math:`p_{ij}^{(2)} = -\lambda_3(D_{ij}-2) + \lambda_4(D_{ij}-2)^2`. The total penalty term for any pair of beads is :math:`p_{ij} = p_{ij}^{(1)} + p_{ij}^{(2)}`.

        Returns:
            A qubit operator for the H_olap term.
        """
        penalty_olap_1 = self._penalty_parameters.penalty_olap_1
        penalty_olap_2 = self._penalty_parameters.penalty_olap_2
        penalty_olap_3 = self._penalty_parameters.penalty_olap_3
        penalty_olap_4 = self._penalty_parameters.penalty_olap_4
        num_qubits = self._distance_map._num_qubits
        full_id = build_full_identity(num_qubits)
        h_olap = 0
        for i in range(self._peptide_length - 4):
            for j in range(
                i + 4, self._peptide_length
            ):  # overlap cannot happen within 3 beads
                distance = self._distance_map.distance_map[i][j]
                h_olap += -penalty_olap_1 * ((2 * (j - i) ** 2) * full_id - distance)
                h_olap += (
                    penalty_olap_2 * ((2 * (j - i) ** 2) * full_id - distance) ** 2
                )
                h_olap += -penalty_olap_3 * (distance - 2 * full_id)
                h_olap += penalty_olap_4 * (distance - 2 * full_id) ** 2
        return fix_qubits(h_olap)

    def _create_h_contact(self) -> SparsePauliOp:
        """
        Creates qubit operators for the H_contact term for first nearest neighbor interactions.

        Returns:
            A qubit operator for the H_contact term.
        """
        h_contact = 0
        for i in range(self._peptide_length - 2):
            for j in range(i + 2, self._peptide_length):
                h_contact += (self._contact_map.contact_map[i][j]) ^ (
                    self._distance_map.first_neighbor(i, j, self._pair_energies)
                )
        return fix_qubits(h_contact)
