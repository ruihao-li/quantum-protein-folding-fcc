# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

"""A class defining the Miyazawa-Jernigan interaction between beads of a peptide."""

import os
import numpy as np


def _load_energy_matrix_file(
    file_name: str = "mj_matrix",
) -> tuple[np.ndarray, list[str]]:
    """Returns the energy matrix from the Miyazawa-Jernigan potential file."""

    path = _construct_resource_path(file_name)
    matrix = np.loadtxt(fname=path, dtype=str)
    energy_matrix = _parse_energy_matrix(matrix)
    symbols = list(matrix[0, :])
    return energy_matrix, symbols


def _construct_resource_path(file_name: str = "mj_matrix") -> str:
    path = os.path.realpath(
        os.path.join(
            os.path.dirname(__file__),
            "potentials",
            file_name + ".txt",
        )
    )

    return os.path.normpath(path)


def _parse_energy_matrix(matrix: np.ndarray) -> np.ndarray:
    """Parses a matrix loaded from the Miyazawa-Jernigan potential file."""
    energy_matrix = np.zeros((np.shape(matrix)[0], np.shape(matrix)[1]))
    for row in range(1, np.shape(matrix)[0]):
        for col in range(row - 1, np.shape(matrix)[1]):
            energy_matrix[row, col] = float(matrix[row, col])
    energy_matrix = energy_matrix[1:,]
    return energy_matrix


def _validate_residue_sequence(residue_sequence: str):
    """
    Checks if the provided residue sequence contains allowed characters.

    Args:
        residue_sequence: A list or a string that contains characters defining
        residues for a chain of proteins.

    Raises:
        ValueError: If an illegal residue character is discovered.
    """
    for residue_symbol in residue_sequence:
        _validate_residue_symbol(residue_symbol)


def _validate_residue_symbol(residue_symbol: str):
    """
    Checks if the provided residue character is legal. If not, an
    ValueError is raised.

    Args:
        residue_symbol: symbol of a residue.

    Raises:
        ValueError: If a symbol provided is not legal.
    """
    valid_residues = [
        "A",  # Alanine
        "C",  # Cysteine
        "D",  # Aspartic acid
        "E",  # Glutamic acid
        "F",  # Phenylalanine
        "G",  # Glycine
        "H",  # Histidine
        "I",  # Isoleucine
        "K",  # Lysine
        "L",  # Leucine
        "M",  # Methionine
        "N",  # Asparagine
        "P",  # Proline
        "Q",  # Glutamine
        "R",  # Arginine
        "S",  # Serine
        "T",  # Threonine
        "V",  # Valine
        "W",  # Tryptophan
        "Y",  # Tyrosine
    ]
    if residue_symbol != "" and residue_symbol not in valid_residues:
        raise ValueError(
            f"Provided residue type {residue_symbol} is not valid. Valid residue types are "
            f"{valid_residues}"
        )


class MiyazawaJerniganInteraction:
    """A class defining a Miyazawa-Jernigan interaction between beads of a
    peptide. Details of this model can be found in Miyazawa, S. and Jernigan, R.
    L. J. Mol. Biol.256, 623–644 (1996), Table 3."""

    def __init__(self, energy_matrix_file: str = "mj_matrix"):
        """
        Args:
            energy_matrix_file: Name of the file containing the
            Miyazawa-Jernigan potential.
        """
        self.energy_matrix_file = energy_matrix_file

    def validate_residue_sequence(self, residue_sequence: str) -> None:
        """Validate that the sequence matches the Miyazawa-Jernigan alphabet."""
        _validate_residue_sequence(residue_sequence)

    def calculate_energy_matrix(self, residue_sequence: str) -> np.ndarray:
        """
        Calculates an energy matrix for a Miyazawa-Jernigan interaction based on
        the Miyazawa-Jernigan potential file.

        Args:
            residue_sequence: A string that contains characters defining
            residues for a chain of proteins.

        Returns:
            Numpy array of pair energies for amino acids.
        """
        chain_len = len(residue_sequence)
        self.validate_residue_sequence(residue_sequence)
        mj_interaction, list_aa = _load_energy_matrix_file(self.energy_matrix_file)
        pair_energies = np.zeros((chain_len, chain_len))
        for i in range(chain_len):
            for j in range(i + 1, chain_len):
                aa_i = list_aa.index(residue_sequence[i])
                aa_j = list_aa.index(residue_sequence[j])
                pair_energies[i, j] = mj_interaction[min(aa_i, aa_j), max(aa_i, aa_j)]
        return pair_energies
