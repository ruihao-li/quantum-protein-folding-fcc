"""A class for the result of the protein folding problem."""

from .fcc_peptide import Peptide
from .fcc_protein_shape import ProteinShapeDecoder, ProteinShapeFileGen
import numpy as np
import matplotlib.pyplot as plt


class ProteinFoldingResult:

    def __init__(
        self, peptide: Peptide, unused_qubits: list[int], solution_bitstring: str
    ):
        """
        Args:
            peptide: The peptide defining the protein subject to the folding problem.
            unused_qubits: The list of indices for qubits in the original problem formulation that were removed during compression.
            solution_bitstring: The compact bitstring representing both the configuration and interaction qubits.
        """
        self._solution_bitstring = solution_bitstring
        self._unused_qubits = unused_qubits
        self._peptide = peptide

        self._protein_shape_decoder = ProteinShapeDecoder(
            peptide=self._peptide,
            solution_bitstring=self._solution_bitstring,
        )

        self._protein_shape_file_gen = ProteinShapeFileGen(
            peptide=self._peptide,
            turn_sequence=self._protein_shape_decoder.turn_sequence,
        )

    @property
    def protein_shape_decoder(self) -> ProteinShapeDecoder:
        """Returns the :class:`ProteinShapeDecoder` of the result.
        This class will interpret the result bitstring and return the encoded information.
        """
        return self._protein_shape_decoder

    @property
    def protein_shape_file_gen(self) -> ProteinShapeFileGen:
        """Returns the :class:`ProteinShapeFileGen` of the result."""
        return self._protein_shape_file_gen

    @property
    def solution_bitstring(self) -> str:
        """Returns the bitstring encoding the turns of the protein."""
        return self._solution_bitstring

    @property
    def turn_sequence(self) -> list:
        """Returns the turns of the protein."""
        return self.protein_shape_decoder.turn_sequence

    def get_result_binary_vector(self) -> str:
        """Returns the bitstring that encodes the solution of a protein folding problem."""
        unused_qubits = self._unused_qubits
        result = []
        offset = 0
        size = len(self._solution_bitstring)
        for i in range(size):
            index = size - 1 - i
            while i + offset in unused_qubits:
                result.append("_")
                offset += 1
            result.append(self._solution_bitstring[index])

        return "".join(result[::-1])

    def save_xyz_file(
        self,
        name: str | None = None,
        path: str = "",
        comment: str = "",
        replace: bool = False,
    ) -> None:
        """
        Generates and saves a .xyz file.

        Args:
            name: Name of the file to be generated. If the name is ``None`` the name of the file will be the letters of the amino acids on the peptide chain. If a file of the same name already exists then the action taken is dependent on the `replace` arg.
            path: Path where the file will be generated. If left empty the file will be saved in the working directory.
            comment: Comment to be added to the second line of the file. By default, the line will be left blank.
            replace: If ``True``, the file will be overwritten if it already exists.

        Raises:
            FileExistsError: If the file already exists and replace is ``False``.
        """
        if name is None:
            name = str(self._peptide.peptide_sequence)
        self.protein_shape_file_gen.save_xyz_file(
            name=name, path=path, comment=comment, replace=replace
        )

    def get_figure(
        self, title: str = "Protein Structure", ticks: bool = False, grid: bool = False
    ) -> plt.Figure:
        """
        Generates a figure of the peptide in 3D.

        Args:
            title: The title of the plot.
            ticks: Boolean for showing ticks in the graphic.
            grid: Boolean for showing the grid in the graphic.
        Returns:
            A figure with the folded protein.
        """
        (
            _x_main,
            _y_main,
            _z_main,
        ) = np.split(self.protein_shape_file_gen.amino_acid_positions.transpose(), 3, 0)
        _x_main, _y_main, _z_main = (
            _x_main[0],
            _y_main[0],
            _z_main[0],
        )
        fig = plt.figure()
        ax_graph = fig.add_subplot(projection="3d")
        for i, amino_acid in enumerate(self.protein_shape_file_gen.amino_acid_list):
            ax_graph.text(
                _x_main[i],
                _y_main[i],
                _z_main[i],
                amino_acid,
                size=10,
                zorder=10,
                color="k",
            )

        ax_graph.plot3D(_x_main, _y_main, _z_main)
        ax_graph.scatter3D(_x_main, _y_main, _z_main, s=500)
        ax_graph.set_box_aspect([1, 1, 1])

        ax_graph.grid(grid)

        if not ticks:
            ax_graph.set_xticks([])
            ax_graph.set_yticks([])
            ax_graph.set_zticks([])

        ax_graph.set_xlabel("x")
        ax_graph.set_ylabel("y")
        ax_graph.set_zlabel("z")

        ax_graph.set_title(title)
        return fig
