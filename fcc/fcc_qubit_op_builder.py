# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

"""Builds the total qubit operator for the full Hamiltonian encoding a protein folding problem on the FCC lattice."""

import numpy as np
from sklearn.metrics import r2_score
from qiskit.quantum_info import SparsePauliOp
from .fcc_peptide import Peptide
from .fcc_distance_map import DistanceMap
from .fcc_contact_map import ContactMap
from .fcc_penalty_parameters import PenaltyParameters
from .utils import (
    build_full_identity,
    fix_qubits,
    compose_IZ_ops,
    remove_unused_qubits,
)
import time


class QubitOpBuilder:

    def __init__(
        self,
        peptide: Peptide,
        pair_energies: np.ndarray,
        penalty_parameters: PenaltyParameters,
        *,
        build_geometry_maps: bool = True,
    ):
        """
        Args:
            peptide: A Peptide object that includes all information about a protein.
            pair_energies: A matrix of pair energies between beads.
            penalty_parameters: Parameters that define the strength of constraints enforcing in the problem.
            build_geometry_maps: Whether to construct the symbolic distance and
                contact maps eagerly. The turn-only path disables this.
        """
        self._peptide = peptide
        self._peptide_length = peptide.peptide_length
        self._pair_energies = pair_energies
        self._penalty_parameters = penalty_parameters
        self._contact_map = ContactMap(peptide) if build_geometry_maps else None
        self._distance_map = DistanceMap(peptide) if build_geometry_maps else None
        self._num_config_qubits = 4 * (self._peptide_length - 1)
        self._num_contact_qubits = (
            self._peptide_length**2 - 3 * self._peptide_length + 2
        ) // 2

    def _ensure_geometry_maps(self) -> None:
        """Construct the legacy symbolic geometry maps on first use."""

        if self._contact_map is None:
            self._contact_map = ContactMap(self._peptide)
        if self._distance_map is None:
            self._distance_map = DistanceMap(self._peptide)

    def build_qubit_op(self, r2_threshold: float, chunk: int = 20) -> SparsePauliOp:
        """
        Builds the total qubit operator for the full Hamiltonian encoding a
        protein folding problem. H_total = H_back + H_redun + H_olap +
        H_contact.

        Args:
            r2_threshold: The threshold for the R^2 score of the Chebyshev fit
            when building the non-overlapping constraint.
            chunk: Size of the chunks to split the qubit operator into.

        Returns:
            A qubit operator for the full Hamiltonian encoding a protein folding
            problem.
        """
        self._ensure_geometry_maps()
        contact_id = build_full_identity(self._num_contact_qubits)
        h_back = self._create_h_back()
        if h_back != 0:
            h_back = contact_id ^ h_back
        h_redun = self._create_h_redun()
        if h_redun != 0:
            h_redun = contact_id ^ h_redun
        h_olap = self._create_h_olap(r2_threshold, chunk)
        if h_olap != 0:
            h_olap = contact_id ^ h_olap
        h_contact = self._create_h_contact()
        h_total = h_back + h_redun + h_olap + h_contact
        return h_total.simplify()

    def build_olap_constr_ops(self) -> dict[tuple[int, int], SparsePauliOp]:
        """
        Builds qubit operators for the constraints that penalize overlapping
        beads, which are subsequently used in the VQEC approach based on the
        Lagrangian dual method (arXiv:2311.08502). Note that this should be used
        only when the overlap penalty parameter is None.

        Returns:
            A dictionary containing the indices of bead pairs as keys and the
            corresponding qubit operators as values.
        """
        self._ensure_geometry_maps()
        contact_id = build_full_identity(self._num_contact_qubits)
        olap_constraints = self._create_olap_constraints()
        for key in olap_constraints:
            olap_constraints[key] = contact_id ^ olap_constraints[key]
        return olap_constraints

    def build_turn_penalty_ops(self) -> tuple[SparsePauliOp, SparsePauliOp]:
        """
        Builds compact backtracking and redundant-code Hamiltonians.

        The six symmetry-fixed turn qubits are removed without constructing
        distance maps, contact maps, interaction ancillas, or overlap fits.

        Returns:
            A tuple containing the compact backtracking Hamiltonian
            followed by the compact redundant-code Hamiltonian.

        Raises:
            ValueError: If the peptide contains fewer than three residues.
        """

        if self._peptide_length < 3:
            raise ValueError(
                "Compact FCC turn Hamiltonians require at least three residues"
            )
        compact_qubits = 4 * self._peptide_length - 10
        fixed_qubits = [0, 1, 2, 3, 6, 7]
        zero = SparsePauliOp(
            "I" * compact_qubits, coeffs=np.asarray([0.0], dtype=complex)
        )

        compact_ops = []
        for operator in (self._create_h_back(), self._create_h_redun()):
            if not isinstance(operator, SparsePauliOp):
                compact_ops.append(zero)
                continue
            compact, _ = remove_unused_qubits(operator, fixed_qubits)
            compact_ops.append(compact.simplify())
        return compact_ops[0], compact_ops[1]

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
        Creates qubit operators for the H_back term, which penalizes consecutive
        turns in opposite directions. Note that the first bead (index 0) is
        omitted because the two turns following it cannot be in opposite
        directions by construction.

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
        Creates qubit operators for the H_redun term, which penalizes redundant
        turns that do not correspond to any physical turns. Note that the first
        two turns are omitted because they are always physical by construction.

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

    def _compose_in_chunks(
        self, op1: SparsePauliOp, op2: SparsePauliOp, chunk: int = 20
    ) -> SparsePauliOp:
        """
        Composes two qubit operators by splitting the second operator into chunks (to avoid memory errors).

        Args:
            op1: The first qubit operator.
            op2: The second qubit operator.
            chunk: Size of the chunks to split the second operator into.

        Returns:
            The composed qubit operator.
        """
        final_op = 0

        for i in range(0, op2.size, chunk):
            terms = op2[i : i + chunk]
            final_op += compose_IZ_ops(op1, terms)
            final_op = final_op.simplify()

        return final_op

    def _create_h_olap(self, r2_threshold: float, chunk: int = 20) -> SparsePauliOp:
        r"""
        Creates qubit operators for the H_olap term, which penalizes overlapping
        beads. To ensure the non-overlapping condition, we impose constraints on
        the distance function that encodes the squared distance between beads,
        that is, for any two beads :math:`i` and :math:`j`, we require :math:`2
        \leq D_{ij} \leq 2(i-j)^2`. For each pair of beads, we define
        :math:`h_{ij} = D_{ij} - 2`. Ideally, the penalty function would have a
        large positive value (set by `penalty_olap`) when :math:`h_{ij} = -2`
        (corresponding to the case when the beads overlap) and zero otherwise.
        Practically, we cannot enforce this constraint directly, so we use
        polynomials to approximate the penalty function. Note that the degree of
        the polynomial required to achieve an accuracy threshold
        :math:`R^2_{ij}` depends on the domain of the function, which is
        :math:`[-2, 2(i-j)^2 - 2]`. The further the pair of beads are from each
        other, the higher the degree of the polynomial required to approximate
        the penalty function. Note that if the overlap penalty parameter is
        None, this returns 0.

        Args:
            r2_threshold: The threshold for the R^2 score of the Chebyshev fit
            when building the non-overlapping constraint.
            chunk: Size of the chunks to split the second operator into.

        Returns:
            A qubit operator for the H_olap term.
        """
        penalty_olap = self._penalty_parameters.penalty_olap
        if penalty_olap is None:
            return 0
        full_id = build_full_identity(self._num_config_qubits)
        h_olap = 0
        for i in range(self._peptide_length - 3):
            for j in range(
                i + 3, self._peptide_length
            ):  # overlap cannot happen within 3 beads on FCC lattice
                dist_op = self._distance_map[(i, j)]
                print(f"Beads {i} and {j} -- dist_op size: {dist_op.size}")
                # Perform Chebyshev fit to approximate the penalty function
                x = np.arange(0, 2 * (j - i) ** 2 + 1, 2)  # x = D_{ij}
                y = [penalty_olap] + [0] * (len(x) - 1)
                # Initialize the max degree of the polynomial
                degree = 2
                cheb_fit = np.polynomial.Chebyshev.fit(x, y, degree)(x)
                r2 = r2_score(y, cheb_fit)
                while r2 < r2_threshold:
                    degree += 1
                    cheb_fit = np.polynomial.Chebyshev.fit(x, y, degree)
                    r2 = r2_score(y, cheb_fit(x))
                # Convert the Chebyshev fit coefficients to polynomial coefficients
                cheb_coeffs = cheb_fit.convert().coef
                # print(f"Highest degree of the polynomial for {i} and {j}: {degree}")
                poly_coeffs = np.polynomial.chebyshev.cheb2poly(cheb_coeffs)
                # Create the qubit operator based on the polynomial coefficients
                h_olap += poly_coeffs[0] * full_id

                # Perform the composition of operators in chunks
                for k in range(1, len(poly_coeffs)):
                    if k == 1:
                        h_op = dist_op.simplify()
                    else:
                        tic = time.time()
                        h_op = self._compose_in_chunks(h_op, dist_op, chunk)
                        toc = time.time()
                        # print(f"Time taken to do composition @ k = {k}: {toc - tic}")
                    # print(f"h_op size @ k = {k}: {h_op.size}")
                    h_olap = SparsePauliOp.sum(
                        [h_olap, poly_coeffs[k] * h_op]
                    ).simplify()
        return fix_qubits(h_olap)

    def _create_h_contact(self) -> SparsePauliOp:
        """
        Creates qubit operators for the H_contact term for first nearest
        neighbor interactions.

        Returns:
            A qubit operator for the H_contact term.
        """
        h_contact = 0
        for i in range(self._peptide_length - 2):
            for j in range(i + 2, self._peptide_length):
                # contact_map contains operators that act only on the
                # interaction qubits; distance_map.first_neighbor returns the
                # qubit operator that acts on the configuration qubits
                h_contact += (self._contact_map.contact_map[i][j]) ^ (
                    self._distance_map.first_neighbor(
                        i, j, self._pair_energies, pair_energies_multiplier=1.0
                    )
                )
        return fix_qubits(h_contact)

    def _create_olap_constraints(self) -> dict[tuple[int, int], SparsePauliOp]:
        r"""
        Creates qubit operators for the constraints that penalize overlapping
        beads, which are subsequently used in the VQEC approach based on
        Lagrange multipliers (arXiv:2311.08502). Note that the constraints are
        put in the form of :math:`F_{ij} \leq 0`, where :math:`F_{ij} = 2 -
        D_{ij}`.

        Returns:
            A dictionary containing the indices of bead pairs as keys and the
            corresponding qubit operators as values.
        """
        olap_constraints = {}
        for i in range(self._peptide_length - 3):
            for j in range(i + 3, self._peptide_length):
                dist_op = self._distance_map[(i, j)]
                olap_constraints[(i, j)] = (
                    2 * build_full_identity(self._num_config_qubits) - dist_op
                ).simplify()
        return fix_qubits(olap_constraints)
