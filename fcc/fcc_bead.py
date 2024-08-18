"""A class defining a bead of a peptide (main chain only)."""

from qiskit.quantum_info import SparsePauliOp
from .utils import build_full_identity


class Bead:

    def __init__(
        self,
        main_index: int,
        residue_type: str,
        turn_qubits: tuple[SparsePauliOp, SparsePauliOp, SparsePauliOp, SparsePauliOp],
    ):
        """
        Args:
            main_index: Index of the bead on the main chain in a peptide.
            residue_type: A character representing the type of a residue for the
            bead.
            turn_qubits: A tuple of four Pauli operators that encodes the turn
            following from a given bead index. Note that four qubits are needed
            to encode each turn in the FCC lattice.
        """
        self._main_index = main_index
        self._residue_type = residue_type
        self._turn_qubits = turn_qubits
        if self._residue_type and self._turn_qubits is not None:
            self._full_id = build_full_identity(turn_qubits[0].num_qubits)
            self._turn_indicator_pxpy = self._build_turn_indicator_pxpy()
            self._turn_indicator_pxmy = self._build_turn_indicator_pxmy()
            self._turn_indicator_mxpy = self._build_turn_indicator_mxpy()
            self._turn_indicator_mxmy = self._build_turn_indicator_mxmy()
            self._turn_indicator_pypz = self._build_turn_indicator_pypz()
            self._turn_indicator_pymz = self._build_turn_indicator_pymz()
            self._turn_indicator_mypz = self._build_turn_indicator_mypz()
            self._turn_indicator_mymz = self._build_turn_indicator_mymz()
            self._turn_indicator_pxpz = self._build_turn_indicator_pxpz()
            self._turn_indicator_pxmz = self._build_turn_indicator_pxmz()
            self._turn_indicator_mxpz = self._build_turn_indicator_mxpz()
            self._turn_indicator_mxmz = self._build_turn_indicator_mxmz()
            self._turn_indicator_0010 = self._build_turn_indicator_0010()
            self._turn_indicator_0001 = self._build_turn_indicator_0001()
            self._turn_indicator_1101 = self._build_turn_indicator_1101()
            self._turn_indicator_1110 = self._build_turn_indicator_1110()

    @property
    def turn_qubits(
        self,
    ) -> tuple[SparsePauliOp, SparsePauliOp, SparsePauliOp, SparsePauliOp]:
        """Returns the list of four qubits that encode the turn following from the bead."""
        return self._turn_qubits

    @property
    def residue_type(self) -> str:
        """Returns a residue type."""
        return self._residue_type

    @property
    def main_index(self) -> int:
        """Returns the main index of the bead."""
        return self._main_index

    def _build_turn_indicator_pxpy(self) -> SparsePauliOp:
        """Turn indicator function for (+x, +y) direction: 0000."""
        return (
            (self._full_id - self._turn_qubits[0])
            @ (self._full_id - self._turn_qubits[1])
            @ (self._full_id - self._turn_qubits[2])
            @ (self._full_id - self._turn_qubits[3])
        ).simplify()

    def _build_turn_indicator_pxmy(self) -> SparsePauliOp:
        """Turn indicator function for (+x, -y) direction: 1111."""
        return (
            self._turn_qubits[0]
            @ self._turn_qubits[1]
            @ self.turn_qubits[2]
            @ self.turn_qubits[3]
        ).simplify()

    def _build_turn_indicator_mxpy(self) -> SparsePauliOp:
        """Turn indicator function for (-x, +y) direction: 1100."""
        return (
            self._turn_qubits[0]
            @ self._turn_qubits[1]
            @ (self._full_id - self._turn_qubits[2])
            @ (self._full_id - self._turn_qubits[3])
        ).simplify()

    def _build_turn_indicator_mxmy(self) -> SparsePauliOp:
        """Turn indicator function for (-x, -y) direction: 0011."""
        return (
            (self._full_id - self._turn_qubits[0])
            @ (self._full_id - self._turn_qubits[1])
            @ self._turn_qubits[2]
            @ self._turn_qubits[3]
        ).simplify()

    def _build_turn_indicator_pypz(self) -> SparsePauliOp:
        """Turn indicator function for (+y, +z) direction: 1001."""
        return (
            self._turn_qubits[0]
            @ (self._full_id - self._turn_qubits[1])
            @ (self._full_id - self._turn_qubits[2])
            @ self._turn_qubits[3]
        ).simplify()

    def _build_turn_indicator_pymz(self) -> SparsePauliOp:
        """Turn indicator function for (+y, -z) direction: 1010."""
        return (
            self._turn_qubits[0]
            @ (self._full_id - self._turn_qubits[1])
            @ self._turn_qubits[2]
            @ (self._full_id - self._turn_qubits[3])
        ).simplify()

    def _build_turn_indicator_mypz(self) -> SparsePauliOp:
        """Turn indicator function for (-y, +z) direction: 0110."""
        return (
            (self._full_id - self._turn_qubits[0])
            @ self._turn_qubits[1]
            @ self._turn_qubits[2]
            @ (self._full_id - self._turn_qubits[3])
        ).simplify()

    def _build_turn_indicator_mymz(self) -> SparsePauliOp:
        """Turn indicator function for (-y, -z) direction: 0101."""
        return (
            (self._full_id - self._turn_qubits[0])
            @ self._turn_qubits[1]
            @ (self._full_id - self._turn_qubits[2])
            @ self._turn_qubits[3]
        ).simplify()

    def _build_turn_indicator_pxpz(self) -> SparsePauliOp:
        """Turn indicator function for (+x, +z) direction: 1000."""
        return (
            self._turn_qubits[0]
            @ (self._full_id - self._turn_qubits[1])
            @ (self._full_id - self._turn_qubits[2])
            @ (self._full_id - self._turn_qubits[3])
        ).simplify()

    def _build_turn_indicator_pxmz(self) -> SparsePauliOp:
        """Turn indicator function for (+x, -z) direction: 1011."""
        return (
            self._turn_qubits[0]
            @ (self._full_id - self._turn_qubits[1])
            @ self._turn_qubits[2]
            @ self._turn_qubits[3]
        ).simplify()

    def _build_turn_indicator_mxpz(self) -> SparsePauliOp:
        """Turn indicator function for (-x, +z) direction: 0111."""
        return (
            (self._full_id - self._turn_qubits[0])
            @ self._turn_qubits[1]
            @ self._turn_qubits[2]
            @ self._turn_qubits[3]
        ).simplify()

    def _build_turn_indicator_mxmz(self) -> SparsePauliOp:
        """Turn indicator function for (-x, -z) direction: 0100."""
        return (
            (self._full_id - self._turn_qubits[0])
            @ self._turn_qubits[1]
            @ (self._full_id - self._turn_qubits[2])
            @ (self._full_id - self._turn_qubits[3])
        ).simplify()

    def _build_turn_indicator_0010(self) -> SparsePauliOp:
        """Turn indicator function for the unused state 0010."""
        return (
            (self._full_id - self._turn_qubits[0])
            @ (self._full_id - self._turn_qubits[1])
            @ self._turn_qubits[2]
            @ (self._full_id - self._turn_qubits[3])
        ).simplify()

    def _build_turn_indicator_0001(self) -> SparsePauliOp:
        """Turn indicator function for the unused state 0001."""
        return (
            (self._full_id - self._turn_qubits[0])
            @ (self._full_id - self._turn_qubits[1])
            @ (self._full_id - self._turn_qubits[2])
            @ self._turn_qubits[3]
        ).simplify()

    def _build_turn_indicator_1101(self) -> SparsePauliOp:
        """Turn indicator function for the unused state 1101."""
        return (
            self._turn_qubits[0]
            @ self._turn_qubits[1]
            @ (self._full_id - self._turn_qubits[2])
            @ self._turn_qubits[3]
        ).simplify()

    def _build_turn_indicator_1110(self) -> SparsePauliOp:
        """Turn indicator function for the unused state 1110."""
        return (
            self._turn_qubits[0]
            @ self._turn_qubits[1]
            @ self._turn_qubits[2]
            @ (self._full_id - self._turn_qubits[3])
        ).simplify()

    @property
    def physical_turn_indicators(
        self,
    ) -> None | tuple[SparsePauliOp, ...]:
        """Returns all physical turn indicator functions for the bead."""
        if self._turn_qubits is None:
            return None
        return (
            self._turn_indicator_pxpy,
            self._turn_indicator_pxmy,
            self._turn_indicator_mxpy,
            self._turn_indicator_mxmy,
            self._turn_indicator_pypz,
            self._turn_indicator_pymz,
            self._turn_indicator_mypz,
            self._turn_indicator_mymz,
            self._turn_indicator_pxpz,
            self._turn_indicator_pxmz,
            self._turn_indicator_mxpz,
            self._turn_indicator_mxmz,
        )

    @property
    def unused_turn_indicators(
        self,
    ) -> None | tuple[SparsePauliOp, ...]:
        """Returns all unused turn indicator functions for the bead."""
        if self._turn_qubits is None:
            return None
        return (
            self._turn_indicator_0010,
            self._turn_indicator_0001,
            self._turn_indicator_1101,
            self._turn_indicator_1110,
        )
