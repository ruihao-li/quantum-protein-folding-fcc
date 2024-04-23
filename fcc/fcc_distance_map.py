"""A class for building a distance map that stores distances between beads in a peptide on the FCC lattice."""

from collections import defaultdict
import numpy as np
from qiskit.quantum_info import SparsePauliOp
from .fcc_peptide import Peptide
from .utils import fix_qubits, build_full_identity


class DistanceMap:

    def __init__(self, peptide: Peptide):
        """
        Args:
            peptide: A Peptide object that includes all information about a protein.
        """
        self._peptide = peptide
        self._peptide_length = peptide.peptide_length
        self._num_qubits = 4 * (self._peptide_length - 1)
        (self._distance_map, self._num_distances) = self._build_distance_map()

    @property
    def peptide(self) -> Peptide:
        """Returns a peptide."""
        return self._peptide

    @property
    def distance_map(self) -> defaultdict[int, dict[int, SparsePauliOp]]:
        """Returns a distance map."""
        return self._distance_map

    @property
    def num_distances(self) -> int:
        """Returns the number of distances calculated."""
        return self._num_distances

    def _build_distance_map(
        self,
    ) -> tuple[defaultdict[int, dict[int, SparsePauliOp]], int]:
        """
        Builds a distance map for a given peptide, which contains the squared distances between all pairs of beads on the main chain.

        Args:
            peptide: A Peptide object that includes all information about a protein.

        Returns:
            A tuple of a distance map and the number of distances calculated.
        """
        num_distances = 0
        distance_map = defaultdict(dict)
        for lower_bead_idx in range(self._peptide_length - 1):
            for upper_bead_idx in range(lower_bead_idx + 1, self._peptide_length):
                distance_map[lower_bead_idx][upper_bead_idx] = self._compute_distance(
                    lower_bead_idx, upper_bead_idx
                )
                num_distances += 1
        return distance_map, num_distances

    def _compute_x_coordinate(self, bead_idx: int) -> SparsePauliOp:
        """
        Computes the x-coordinate of a bead on the FCC lattice.

        Args:
            bead_idx: Index of the bead in the peptide.

        Returns:
            A Pauli operator that encodes the x-coordinate of the bead.

        Raises:
            ValueError: If the bead index is out of range.
        """
        if bead_idx >= self._peptide_length:
            raise ValueError(
                f"Bead index cannot be great than {self._peptide_length - 1}."
            )
        if bead_idx == 0:
            return SparsePauliOp.from_list([("I" * self._num_qubits, 0)])
        x_coordinate = 0
        for j in range(0, bead_idx):
            bead = self._peptide.beads_list[j]
            turn_indicator_pxpy = bead._turn_indicator_pxpy
            turn_indicator_pxmy = bead._turn_indicator_pxmy
            turn_indicator_pxpz = bead._turn_indicator_pxpz
            turn_indicator_pxmz = bead._turn_indicator_pxmz
            turn_indicator_mxpy = bead._turn_indicator_mxpy
            turn_indicator_mxmy = bead._turn_indicator_mxmy
            turn_indicator_mxpz = bead._turn_indicator_mxpz
            turn_indicator_mxmz = bead._turn_indicator_mxmz
            x_coordinate += fix_qubits(
                turn_indicator_pxpy
                + turn_indicator_pxmy
                + turn_indicator_pxpz
                + turn_indicator_pxmz
                - turn_indicator_mxpy
                - turn_indicator_mxmy
                - turn_indicator_mxpz
                - turn_indicator_mxmz
            )
        return x_coordinate.simplify()

    def _compute_y_coordinate(self, bead_idx: int) -> SparsePauliOp:
        """
        Computes the y-coordinate of a bead on the FCC lattice.

        Args:
            bead_idx: Index of the bead in the peptide.

        Returns:
            A Pauli operator that encodes the y-coordinate of the bead.

        Raises:
            ValueError: If the bead index is out of range.
        """
        if bead_idx >= self._peptide_length:
            raise ValueError(
                f"Bead index cannot be great than {self._peptide_length - 1}."
            )
        if bead_idx == 0:
            return SparsePauliOp.from_list([("I" * self._num_qubits, 0)])
        y_coordinate = 0
        for j in range(0, bead_idx):
            bead = self._peptide.beads_list[j]
            turn_indicator_pxpy = bead._turn_indicator_pxpy
            turn_indicator_mxpy = bead._turn_indicator_mxpy
            turn_indicator_pypz = bead._turn_indicator_pypz
            turn_indicator_pymz = bead._turn_indicator_pymz
            turn_indicator_pxmy = bead._turn_indicator_pxmy
            turn_indicator_mxmy = bead._turn_indicator_mxmy
            turn_indicator_mypz = bead._turn_indicator_mypz
            turn_indicator_mymz = bead._turn_indicator_mymz
            y_coordinate += fix_qubits(
                turn_indicator_pxpy
                + turn_indicator_mxpy
                + turn_indicator_pypz
                + turn_indicator_pymz
                - turn_indicator_pxmy
                - turn_indicator_mxmy
                - turn_indicator_mypz
                - turn_indicator_mymz
            )
        return y_coordinate.simplify()

    def _compute_z_coordinate(self, bead_idx: int) -> SparsePauliOp:
        """
        Computes the z-coordinate of a bead on the FCC lattice.

        Args:
            bead_idx: Index of the bead in the peptide.

        Returns:
            A Pauli operator that encodes the z-coordinate of the bead.

        Raises:
            ValueError: If the bead index is out of range.
        """
        if bead_idx >= self._peptide_length:
            raise ValueError(
                f"Bead index cannot be great than {self._peptide_length - 1}."
            )
        if bead_idx == 0:
            return SparsePauliOp.from_list([("I" * self._num_qubits, 0)])
        z_coordinate = 0
        for j in range(0, bead_idx):
            bead = self._peptide.beads_list[j]
            turn_indicator_pxpz = bead._turn_indicator_pxpz
            turn_indicator_mxpz = bead._turn_indicator_mxpz
            turn_indicator_pypz = bead._turn_indicator_pypz
            turn_indicator_mypz = bead._turn_indicator_mypz
            turn_indicator_pxmz = bead._turn_indicator_pxmz
            turn_indicator_mxmz = bead._turn_indicator_mxmz
            turn_indicator_pymz = bead._turn_indicator_pymz
            turn_indicator_mymz = bead._turn_indicator_mymz
            z_coordinate += fix_qubits(
                turn_indicator_pxpz
                + turn_indicator_mxpz
                + turn_indicator_pypz
                + turn_indicator_mypz
                - turn_indicator_pxmz
                - turn_indicator_mxmz
                - turn_indicator_pymz
                - turn_indicator_mymz
            )
        return z_coordinate.simplify()

    def _compute_distance(
        self, lower_bead_idx: int, upper_bead_idx: int
    ) -> SparsePauliOp:
        """
        Computes the distance squared between any two beads on the FCC lattice.

        Args:
            lower_bead_idx: Index of the lower bead in the peptide (assumed to be smaller than the upper bead index).
            upper_bead_idx: Index of the upper bead in the peptide (assumed to be greater than the lower bead index).

        Returns:
            A Pauli operator that encodes the distance squared between two beads.
        """
        x_lower = self._compute_x_coordinate(lower_bead_idx)
        y_lower = self._compute_y_coordinate(lower_bead_idx)
        z_lower = self._compute_z_coordinate(lower_bead_idx)
        x_upper = self._compute_x_coordinate(upper_bead_idx)
        y_upper = self._compute_y_coordinate(upper_bead_idx)
        z_upper = self._compute_z_coordinate(upper_bead_idx)
        return fix_qubits(
            (x_upper - x_lower) ** 2
            + (y_upper - y_lower) ** 2
            + (z_upper - z_lower) ** 2
        )

    def first_neighbor(
        self,
        lower_bead_idx: int,
        upper_bead_idx: int,
        pair_energies: np.ndarray,
        pair_energies_multiplier: float = 1,
    ) -> SparsePauliOp:
        """
        Creates the first nearest neighbor interaction between two beads on the FCC lattice if they are at squared distance of 2 units from each other.

        Args:
            lower_bead_idx: Index of the lower bead in the peptide.
            upper_bead_idx: Index of the upper bead in the peptide.
            pair_energies: Numpy array of pair energies for amino acids; these pair energies are negative.
            pair_energies_multiplier: A constant that multiplies pair energy contributions.

        Returns:
            Contribution to an energetic Hamiltonian (without interaction qubits).
        """
        energy = pair_energies[lower_bead_idx][upper_bead_idx]
        distance_op = self.distance_map[lower_bead_idx][upper_bead_idx]
        # (3 - D) * E; D = 0 is penalized by the non-overlapping constraint; D > 2 has higher energies because E < 0
        expression = (
            pair_energies_multiplier
            * energy
            * (3 * build_full_identity(distance_op.num_qubits) - distance_op)
        )
        return fix_qubits(expression)
