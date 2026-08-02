# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

"""FCC turn encoding, decoding, and coordinate-generation utilities."""

import os
from collections.abc import Iterable, Sequence

import numpy as np

from .fcc_peptide import Peptide


FCC_TURN_VECTORS = np.asarray(
    [
        (1, 1, 0),
        (-1, -1, 0),
        (-1, 1, 0),
        (1, -1, 0),
        (0, 1, 1),
        (0, -1, -1),
        (0, 1, -1),
        (0, -1, 1),
        (1, 0, 1),
        (-1, 0, -1),
        (1, 0, -1),
        (-1, 0, 1),
    ],
    dtype=np.int8,
)
FCC_TURN_VECTORS.setflags(write=False)

TURN_BITS_TO_INDEX = {
    "0000": 0,
    "0011": 1,
    "1100": 2,
    "1111": 3,
    "1001": 4,
    "0101": 5,
    "1010": 6,
    "0110": 7,
    "1000": 8,
    "0100": 9,
    "1011": 10,
    "0111": 11,
    "0001": None,
    "0010": None,
    "1101": None,
    "1110": None,
}
TURN_INDEX_TO_BITS = {
    turn: bits for bits, turn in TURN_BITS_TO_INDEX.items() if turn is not None
}
TURN_INDEX_TO_LABEL = tuple("0123456789ab")


def compact_turn_qubit_count(peptide_length: int) -> int:
    """Return the compact FCC register size ``4N - 10``."""

    length = int(peptide_length)
    if length < 3:
        raise ValueError("The compact FCC encoding requires at least three residues")
    return 4 * length - 10


def compact_turn_qubit_blocks(peptide_length: int) -> tuple[tuple[int, ...], ...]:
    """Return the compact qubit blocks associated with unfixed turn codes."""

    num_qubits = compact_turn_qubit_count(peptide_length)
    blocks = [tuple(range(0, min(2, num_qubits)))]
    blocks.extend(
        tuple(range(start, min(start + 4, num_qubits)))
        for start in range(2, num_qubits, 4)
    )
    return tuple(block for block in blocks if block)


def bitstring_to_turn_sequence(bitstring: str) -> tuple[int | None, ...]:
    """Decode a full turn bitstring in Qiskit's displayed bit order."""

    if len(bitstring) % 4:
        raise ValueError("A full turn bitstring must contain four bits per turn")
    if set(bitstring) - {"0", "1"}:
        raise ValueError("A turn bitstring may contain only '0' and '1'")
    little_endian = bitstring[::-1]
    return tuple(
        TURN_BITS_TO_INDEX[little_endian[offset : offset + 4]]
        for offset in range(0, len(little_endian), 4)
    )


def decode_compact_turn_bitstring(
    compact_bitstring: str, peptide_length: int
) -> tuple[int | None, ...]:
    """Decode a compact configuration bitstring into physical turn order."""

    expected = compact_turn_qubit_count(peptide_length)
    if len(compact_bitstring) != expected:
        raise ValueError(
            f"Expected {expected} compact turn bits, received {len(compact_bitstring)}"
        )
    if set(compact_bitstring) - {"0", "1"}:
        raise ValueError("A compact turn bitstring may contain only '0' and '1'")

    full_bitstring = compact_bitstring + "0000"
    full_bitstring = full_bitstring[:-6] + "00" + full_bitstring[-6:]
    return bitstring_to_turn_sequence(full_bitstring)


def decode_compact_turn_index(
    index: int, peptide_length: int
) -> tuple[int | None, ...]:
    """Decode a compact-register basis index into physical turn order."""

    num_qubits = compact_turn_qubit_count(peptide_length)
    configuration_index = int(index)
    if not 0 <= configuration_index < 2**num_qubits:
        raise ValueError(
            f"Configuration index must lie in [0, {2**num_qubits})"
        )
    return decode_compact_turn_bitstring(
        format(configuration_index, f"0{num_qubits}b"), peptide_length
    )


def format_turn_sequence(
    turn_sequence: Iterable[int | None], *, reverse: bool = False
) -> str:
    """Format turn indices with hexadecimal labels and ``x`` for unused codes."""

    turns = tuple(turn_sequence)
    if reverse:
        turns = turns[::-1]
    labels = []
    for turn in turns:
        if turn is None or int(turn) == -1:
            labels.append("x")
        elif 0 <= int(turn) < len(TURN_INDEX_TO_LABEL):
            labels.append(TURN_INDEX_TO_LABEL[int(turn)])
        else:
            raise ValueError(f"Unknown FCC turn index {turn!r}")
    return "".join(labels)


def turn_sequence_to_lattice_positions(
    turn_sequence: Sequence[int | None], *, invalid_turns_as_zero: bool = False
) -> np.ndarray:
    """Return unscaled integer FCC coordinates for a turn sequence."""

    positions = np.zeros((len(turn_sequence) + 1, 3), dtype=np.int16)
    for index, turn in enumerate(turn_sequence):
        if turn is None or int(turn) == -1:
            if not invalid_turns_as_zero:
                raise ValueError("The turn sequence contains an unused FCC code")
            displacement = np.zeros(3, dtype=np.int16)
        elif 0 <= int(turn) < len(FCC_TURN_VECTORS):
            displacement = FCC_TURN_VECTORS[int(turn)].astype(np.int16)
        else:
            raise ValueError(f"Unknown FCC turn index {turn!r}")
        positions[index + 1] = positions[index] + displacement
    return positions


def noncovalent_residue_pairs(
    peptide_length: int, *, minimum_separation: int = 2
) -> tuple[tuple[int, int], ...]:
    """Return ordered residue pairs separated by at least the given amount."""

    length = int(peptide_length)
    separation = int(minimum_separation)
    if length < 1:
        raise ValueError("peptide_length must be positive")
    if separation < 1:
        raise ValueError("minimum_separation must be positive")
    return tuple(
        (lower, upper)
        for lower in range(length)
        for upper in range(lower + separation, length)
    )


def _compact_bitstring_from_turn_sequence(turn_sequence: Sequence[int]) -> str:
    """Encode turns without validating compact-register symmetry choices."""

    try:
        full_bitstring = "".join(TURN_INDEX_TO_BITS[int(turn)] for turn in turn_sequence)
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Turn sequence contains an invalid FCC turn index") from exc
    compact_bitstring = full_bitstring[:6] + full_bitstring[8:]
    compact_bitstring = compact_bitstring[4:]
    return compact_bitstring[::-1]


def turn_sequence_to_compact_bitstring(turn_sequence: Sequence[int]) -> str:
    """Encode a symmetry-compatible physical turn sequence compactly."""

    turns = tuple(int(turn) for turn in turn_sequence)
    compact_bitstring = _compact_bitstring_from_turn_sequence(turns)
    decoded = decode_compact_turn_bitstring(compact_bitstring, len(turns) + 1)
    if decoded != turns:
        raise ValueError(
            "Turn sequence is incompatible with the compact encoding's fixed symmetry bits"
        )
    return compact_bitstring


def turn_sequence_to_compact_index(turn_sequence: Sequence[int]) -> int:
    """Return the basis index for a symmetry-compatible turn sequence."""

    return int(turn_sequence_to_compact_bitstring(turn_sequence), 2)


class ProteinShapeDecoder:
    """
    This class handles the decoding of the compact solution (bitstring) into the
    turns of the main chain of a protein on the FCC lattice.
    """

    def __init__(self, peptide: Peptide, solution_bitstring: str):
        """
        Args:
            solution_bitstring: The compact bitstring representing both the configuration and interaction qubits.
        """
        self._peptide = peptide
        self._peptide_length = peptide.peptide_length
        self._solution_bitstring = solution_bitstring
        self._turn_sequence = self._get_turn_sequence()
        self._correct_bitstring = self._get_correct_bitstring()

    @property
    def turn_sequence(self) -> list:
        """Returns the turns of the main chain."""
        return self._turn_sequence

    @property
    def correct_bitstring(self) -> str:
        """Returns the correct bitstring for the given configuration qubits."""
        return self._correct_bitstring

    def _bitstring_to_turns(self, bitstring: str) -> list:
        """
        Takes a bitstring encoding the turns of a chain and returns the turns as
        a list of integers. Turn indexing is:
        {
            0: (1, 1, 0)
            1: (-1, -1, 0)
            2: (-1, 1, 0)
            3: (1, -1, 0)
            4: (0, 1, 1)
            5: (0, -1, -1)
            6: (0, 1, -1)
            7: (0, -1, 1)
            8: (1, 0, 1)
            9: (-1, 0, -1)
            10: (1, 0, -1)
            11: (-1, 0, 1)
        }

        Args:
            bitstring: string containing the encoded shape information.
        Returns:
            A list of integers decoding the bitstring.
        """
        return list(bitstring_to_turn_sequence(bitstring))

    def _get_turn_sequence(self) -> list:
        """Returns the turns of the main chain.

        Returns:
            A list of integers representing the sequence of turns of the
            peptide.
        """
        num_qubits = compact_turn_qubit_count(self._peptide_length)
        compact_bitstring = self._solution_bitstring[-num_qubits:]
        return list(
            decode_compact_turn_bitstring(compact_bitstring, self._peptide_length)
        )

    def _get_correct_bitstring(self) -> str:
        """
        The solution bitstring may contain interaction qubit values that are not compatible with the configuration qubits. This returns the correct full bitstring for the given configuration qubits.

        Returns:
            A string representing the correct bitstring for the given
            configuration qubits.
        """
        turn_seq = self._turn_sequence
        # Get the coordinates of the main chain
        coordinates = ProteinShapeFileGen(
            self._peptide, turn_seq
        ).generate_amino_acid_positions()
        # Loop through the pairs of beads and assign the interaction qubits (1 when the beads are in contact and 0 otherwise)
        int_qubits = ""
        for i in range(self._peptide_length - 2):
            for j in range(i + 2, self._peptide_length):
                if np.isclose(np.linalg.norm(coordinates[i] - coordinates[j]), 3.8):
                    int_qubits += "1"
                else:
                    int_qubits += "0"
        int_qubits = int_qubits[::-1]
        config_bitstring = self._solution_bitstring[-(4 * self._peptide_length - 10) :]
        return int_qubits + config_bitstring


class ProteinShapeFileGen:
    """
    This class generates the Cartesian coordinates of the main chain of a
    protein given the turn sequence.
    """

    COORDINATES = 3.8 * (1 / np.sqrt(2)) * FCC_TURN_VECTORS

    def __init__(self, peptide: Peptide, turn_sequence: list):
        """
        Args:
            peptide: The peptide we are getting the positions for.
            turn_sequence: A list of integers encoding the turns of the main
            chain.

        Raises:
            ValueError: If the turn sequence contains None values.
            ValueError: If the length of the turn sequence does not match the
            length of the peptide.
        """
        self._peptide = peptide
        self._peptide_length = peptide.peptide_length
        self._turn_sequence = turn_sequence
        if None in self._turn_sequence:
            raise ValueError(
                "This turn sequence contains one or more turns that are not one of the 12 allowed turns on the FCC lattice. Please rerun the optimization."
            )
        # Check length of turn sequence
        if len(self._turn_sequence) != self._peptide_length - 1:
            raise ValueError(
                "The length of the turn sequence does not match the length of the peptide."
            )
        self._amino_acid_list = np.array(list(peptide.peptide_sequence))
        self._amino_acid_positions = self.generate_amino_acid_positions()

    @property
    def amino_acid_positions(self) -> list:
        """Returns the amino acid positions."""
        return self._amino_acid_positions

    @property
    def amino_acid_list(self) -> np.ndarray:
        """Returns the amino acid list."""
        return self._amino_acid_list

    def generate_amino_acid_positions(self) -> np.ndarray:
        """
        Generates the positions of the amino acids in the main chain.

        Returns:
            A list of arrays with the cartesian coordinates of the amino acids.
        """
        lattice_positions = turn_sequence_to_lattice_positions(self._turn_sequence)
        return 3.8 * (1 / np.sqrt(2)) * lattice_positions

    def get_xyz_data(self) -> np.ndarray:
        """
        Returns the xyz data for the amino acids in the main chain.

        Returns:
            An array containing the amino acid letters and the corresponding xyz
            coordinates. The first column contains the amino acid letters, and
            the next three columns contain the x, y, and z coordinates of the
            amino acids in the main chain.
        """
        xyz_data = np.column_stack([self.amino_acid_list, self.amino_acid_positions])
        return xyz_data

    def save_xyz_file(
        self, filename: str, path: str = "", comment: str = "", replace: bool = False
    ) -> None:
        """
        Saves the data as an .xyz file.

        Args:
            filename: The name of the file to save the data to.
            path: The path to save the file to.
            comment: A comment to add to the second line of the file. By
            default, the line will be left blank.
            replace: Whether to replace the file if it already exists.

        Raises:
            ValueError: If the file already exists and replace is False.
        """
        file_path = os.path.join(path, filename + ".xyz")
        if not replace and os.path.exists(file_path):
            raise ValueError(f"File {file_path} already exists.")
        data = self.get_xyz_data()
        number_of_particles = data.shape[0]
        header = f"{number_of_particles} \n {comment}"
        np.savetxt(
            fname=file_path,
            header=header,
            X=data,
            delimiter=" ",
            fmt="%s",
            comments="",
        )


def turns_to_bitstring(peptide: Peptide, turn_sequence: list) -> str:
    """
    Returns the full bitstring for the given turn sequence.

    Args:
        peptide: The peptide we are getting the positions for.
        turn_sequence: A list of integers representing the sequence of turns of the peptide.

    Returns:
        A string representing the full bitstring for the given turn sequence.
    """
    # Get the coordinates of the main chain
    coordinates = ProteinShapeFileGen(
        peptide, turn_sequence
    ).generate_amino_acid_positions()
    # Loop through the pairs of beads and assign the interaction qubits (1 when the beads are in contact and 0 otherwise)
    peptide_length = len(turn_sequence) + 1
    int_qubits = ""
    for i in range(peptide_length - 2):
        for j in range(i + 2, peptide_length):
            if np.isclose(np.linalg.norm(coordinates[i] - coordinates[j]), 3.8):
                int_qubits += "1"
            else:
                int_qubits += "0"
    int_qubits = int_qubits[::-1]

    config_bitstring = _compact_bitstring_from_turn_sequence(turn_sequence)
    full_bitstring = int_qubits + config_bitstring
    return full_bitstring
