"""Auxiliary classes that translate the bitstrings to turn sequences and the corresponding coordinates."""

import numpy as np
import os
from .fcc_peptide import Peptide


class ProteinShapeDecoder:
    """
    This class handles the decoding of the compact solution (bitstring) into the turns of the main chain of a protein on the FCC lattice.
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

    @property
    def turn_sequence(self) -> list:
        """Returns the turns of the main chain."""
        return self._turn_sequence

    def _bitstring_to_turns(self, bitstring: str) -> list:
        """
        Takes a bitstring encoding the turns of a chain and returns the turns as a list of integers. Turn indexing is:
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
            A list of integers representing the sequence of turns of the peptide.
        """
        # Construct the full turn bitstring
        # Add the first 4 bits corresponding to the fixed first turn (0000)
        full_bitstring = self._solution_bitstring + "0000"
        # Add the two qubits at positions 6 & 7 from the right (00)
        full_bitstring = full_bitstring[:-6] + "00" + full_bitstring[-6:]
        # Decode the bitstring
        return self._bitstring_to_turns(full_bitstring)


class ProteinShapeFileGen:
    """
    This class generates the Cartesian coordinates of the main chain of a protein given the turn sequence.
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
            turn_sequence: A list of integers encoding the turns of the main chain.

        Raises:
            ValueError: If the turn sequence contains None values.
            ValueError: If the length of the turn sequence does not match the length of the peptide.
        """
        self._peptide = peptide
        self._peptide_length = peptide.peptide_length
        self._turn_sequence = turn_sequence
        if None in self._turn_sequence:
            raise ValueError(
                "This is an unphysical solution. Please rerun the optimization."
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

    def generate_amino_acid_positions(self) -> list:
        """
        Generates the positions of the amino acids in the main chain.

        Returns:
            A list of arrays with the cartesian coordinates of the amino acids.
        """
        num_turns = len(self._turn_sequence)
        positions = np.zeros((num_turns + 1, 3), dtype=float)
        for i in range(num_turns):
            positions[i + 1] = positions[i] + self.COORDINATES[self._turn_sequence[i]]
        return positions

    def get_xyz_data(self) -> str:
        """
        Returns the xyz data for the amino acids in the main chain.

        Returns:
            A string with the xyz data for the amino acids in the main chain.
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
            comment: A comment to add to the second line of the file. By default, the line will be left blank.
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
