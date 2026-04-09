# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

"""Auxiliary classes that translate the bitstrings to turn sequences and the corresponding coordinates."""

import numpy as np
import os
from .fcc_peptide import Peptide


FIXED_CONFIG_QUBITS = {0, 1, 2, 3, 6, 7}


def _validate_turn_sequence_values(turn_sequence: list) -> None:
    """Validate that a turn sequence uses FCC turn indices in the range 0..11."""
    invalid_turns = sorted(
        {
            repr(turn)
            for turn in turn_sequence
            if not isinstance(turn, (int, np.integer)) or turn < 0 or turn > 11
        }
    )
    if invalid_turns:
        raise ValueError(
            "turn_sequence contains invalid FCC turn indices: "
            f"[{', '.join(invalid_turns)}]. Allowed values are integers from 0 to 11."
        )


class ProteinShapeDecoder:
    """
    This class handles the decoding of the compact solution (bitstring) into the
    turns of the main chain of a protein on the FCC lattice.
    """

    def __init__(
        self,
        peptide: Peptide,
        solution_bitstring: str,
        unused_qubits: list[int] | None = None,
    ):
        """
        Args:
            solution_bitstring: The measured bitstring representing the
                interaction and configuration qubits. This may be either the
                gauge-fixed full bitstring or the compressed bitstring obtained
                after removing unused qubits from the Hamiltonian.
            unused_qubits: Indices of unused qubits removed during Hamiltonian
                compression. These are used to reconstruct the gauge-fixed full
                bitstring before decoding turns.
        """
        self._peptide = peptide
        self._peptide_length = peptide.peptide_length
        self._solution_bitstring = self._normalize_solution_bitstring(
            solution_bitstring, unused_qubits or []
        )
        self._turn_sequence = self._get_turn_sequence()
        self._correct_bitstring = self._get_correct_bitstring()

    @property
    def turn_sequence(self) -> list:
        """Returns the turns of the main chain."""
        return self._turn_sequence

    @property
    def correct_bitstring(self) -> str:
        """Returns the normalized interaction+configuration bitstring.

        The returned string uses the decoder's gauge-fixed configuration-qubit
        convention rather than reinserting the six fixed FCC gauge qubits.
        """
        return self._correct_bitstring

    def _normalize_solution_bitstring(
        self, solution_bitstring: str, unused_qubits: list[int]
    ) -> str:
        """Return the gauge-fixed bitstring expected by the decoder.

        Measurement bitstrings may include register separators (for example,
        spaces inserted by Qiskit when a circuit uses multiple classical
        registers). Strip those separators before validating and decoding.
        """
        normalized_bitstring = solution_bitstring.replace(" ", "").replace("_", "")
        if set(normalized_bitstring) - {"0", "1"}:
            raise ValueError(
                "solution_bitstring contains non-binary characters: "
                f"{solution_bitstring!r}."
            )

        solution_bitstring = normalized_bitstring
        num_contact_qubits = (
            self._peptide_length**2 - 3 * self._peptide_length + 2
        ) // 2
        num_config_qubits = 4 * (self._peptide_length - 1)
        total_qubits = num_contact_qubits + num_config_qubits
        gauge_fixed_total_qubits = total_qubits - len(FIXED_CONFIG_QUBITS)

        expected_compressed_length = total_qubits - len(unused_qubits)

        if unused_qubits:
            if len(solution_bitstring) == total_qubits:
                full_bitstring = solution_bitstring
            elif len(solution_bitstring) == expected_compressed_length:
                unused_qubit_set = set(unused_qubits)
                reconstructed_little_endian: list[str] = []
                compressed_little_endian = solution_bitstring[::-1]
                compressed_index = 0
                for qubit_index in range(total_qubits):
                    if qubit_index in unused_qubit_set:
                        reconstructed_little_endian.append("0")
                    else:
                        reconstructed_little_endian.append(
                            compressed_little_endian[compressed_index]
                        )
                        compressed_index += 1
                full_bitstring = "".join(reconstructed_little_endian[::-1])
            elif len(solution_bitstring) == gauge_fixed_total_qubits:
                return solution_bitstring
            else:
                raise ValueError(
                    "The provided solution bitstring length is incompatible with the peptide length and unused-qubit map."
                )
        else:
            if len(solution_bitstring) == gauge_fixed_total_qubits:
                return solution_bitstring
            if len(solution_bitstring) == total_qubits:
                full_bitstring = solution_bitstring
            else:
                raise ValueError(
                    "The provided solution bitstring length is incompatible with the peptide length and unused-qubit map."
                )

        full_little_endian = full_bitstring[::-1]
        gauge_fixed_little_endian = [
            bit
            for qubit_index, bit in enumerate(full_little_endian)
            if qubit_index not in FIXED_CONFIG_QUBITS
        ]
        return "".join(gauge_fixed_little_endian[::-1])

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
        # Reverse the bitstring to read from right to left such that the order is configuration qubits followed by the interaction qubits
        bitstring = bitstring[::-1]
        # Unused bitstrings return None
        encoding = {
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
        # Only convert the configuration qubits
        length_turns = len(bitstring) // 4
        return [encoding[bitstring[4 * i : 4 * (i + 1)]] for i in range(length_turns)]

    def _get_turn_sequence(self) -> list:
        """Returns the turns of the main chain.

        Returns:
            A list of integers representing the sequence of turns of the
            peptide.
        """
        # Split the solution bitstring into configuration (4N - 10) and interaction qubits
        config_bitstring = self._solution_bitstring[-(4 * self._peptide_length - 10) :]
        # Add the first 4 bits corresponding to the fixed first turn (0000)
        full_bitstring = config_bitstring + "0000"
        # Add the two qubits at positions 6 & 7 from the right (00)
        full_bitstring = full_bitstring[:-6] + "00" + full_bitstring[-6:]
        # Decode the bitstring
        return self._bitstring_to_turns(full_bitstring)

    def _get_correct_bitstring(self) -> str:
        """
        The solution bitstring may contain interaction-qubit values that are not
        compatible with the configuration qubits. This reconstructs the
        interaction portion while preserving the decoder's normalized
        gauge-fixed configuration-bitstring convention.

        Returns:
            A string representing the normalized bitstring for the given
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

    COORDINATES = (
        3.8
        * (1 / np.sqrt(2))
        * np.array(
            [
                [1, 1, 0],
                [-1, -1, 0],
                [-1, 1, 0],
                [1, -1, 0],
                [0, 1, 1],
                [0, -1, -1],
                [0, 1, -1],
                [0, -1, 1],
                [1, 0, 1],
                [-1, 0, -1],
                [1, 0, -1],
                [-1, 0, 1],
            ]
        )
    )

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
        _validate_turn_sequence_values(self._turn_sequence)
        # Check length of turn sequence
        if len(self._turn_sequence) != self._peptide_length - 1:
            raise ValueError(
                "The length of the turn sequence does not match the length of the peptide."
            )
        self._amino_acid_list = np.array(list(peptide.peptide_sequence))
        self._amino_acid_positions = self.generate_amino_acid_positions()

    @property
    def amino_acid_positions(self) -> np.ndarray:
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
            A NumPy array with the Cartesian coordinates of the amino acids.
        """
        num_turns = len(self._turn_sequence)
        positions = np.zeros((num_turns + 1, 3), dtype=float)
        for i in range(num_turns):
            positions[i + 1] = positions[i] + self.COORDINATES[self._turn_sequence[i]]
        return positions

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
            FileExistsError: If the file already exists and replace is False.
        """
        file_path = os.path.join(path, filename + ".xyz")
        if not replace and os.path.exists(file_path):
            raise FileExistsError(f"File {file_path} already exists.")
        data = self.get_xyz_data()
        number_of_particles = data.shape[0]
        header = f"{number_of_particles}\n{comment}"
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
    Return the gauge-fixed bitstring for a turn sequence.

    The FCC encoding removes four fixed qubits for the first turn and two fixed
    qubits from the second turn. Consequently, only turn sequences whose first
    turn is ``0`` and whose second turn is one of ``{0, 2, 8, 9}`` can be
    represented exactly.

    Args:
        peptide: The peptide we are getting the positions for.
        turn_sequence: A list of integers representing the sequence of turns of the peptide.

    Returns:
        A string representing the gauge-fixed bitstring for the given turn
        sequence.

    Raises:
        ValueError: If the turn sequence is incompatible with the gauge-fixed
            encoding used by the quantum model.
    """
    if len(turn_sequence) != peptide.peptide_length - 1:
        raise ValueError(
            "The length of the turn sequence does not match the length of the peptide."
        )
    if turn_sequence and turn_sequence[0] != 0:
        raise ValueError(
            "The first turn must be 0 because the first FCC turn is fixed by the encoding."
        )
    if len(turn_sequence) > 1 and turn_sequence[1] not in {0, 2, 8, 9}:
        raise ValueError(
            "The second turn must be one of {0, 2, 8, 9} because two of its qubits are fixed by the encoding."
        )

    _validate_turn_sequence_values(turn_sequence)

    turns_bitstring_mapping = {
        0: "0000",
        1: "0011",
        2: "1100",
        3: "1111",
        4: "1001",
        5: "0101",
        6: "1010",
        7: "0110",
        8: "1000",
        9: "0100",
        10: "1011",
        11: "0111",
    }
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

    config_bitstring = "".join(turns_bitstring_mapping[turn] for turn in turn_sequence)
    config_bitstring = config_bitstring[:6] + config_bitstring[8:]
    config_bitstring = config_bitstring[4:]
    config_bitstring = config_bitstring[::-1]
    full_bitstring = int_qubits + config_bitstring
    return full_bitstring
