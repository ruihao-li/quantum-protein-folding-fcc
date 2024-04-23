"""A class for the protein folding problem on an FCC lattice."""

from qiskit.quantum_info import SparsePauliOp
from qiskit_algorithms.minimum_eigensolvers import MinimumEigensolverResult
from fcc_peptide import Peptide
from fcc_penalty_parameters import PenaltyParameters
from fcc_mj_interaction import MiyazawaJerniganInteraction
from fcc_qubit_op_builder import QubitOpBuilder
from fcc_protein_folding_result import ProteinFoldingResult
from utils import remove_unused_qubits


class ProteinFoldingProblem:

    def __init__(
        self,
        peptide: Peptide,
        interaction: MiyazawaJerniganInteraction,
        penalty_parameters: PenaltyParameters,
    ):
        """
        Args:
            peptide: A Peptide object that includes all information about a protein.
            interaction: A Miyazawa-Jernigan interaction object that defines the energy matrix.
            penalty_parameters: Parameters that define the strength of constraints enforcing in the problem.
        """
        self._peptide = peptide
        self._interaction = interaction
        self._penalty_parameters = penalty_parameters
        self._pair_energies = interaction.calculate_energy_matrix(
            peptide.peptide_sequence
        )
        self._qubit_op_builder = QubitOpBuilder(
            self._peptide, self._pair_energies, self._penalty_parameters
        )
        self._unused_qubits = []

    def qubit_op(self) -> SparsePauliOp:
        """
        Builds the total qubit operator for the full Hamiltonian encoding a protein folding problem on the FCC lattice.

        Returns:
            A qubit operator for the full Hamiltonian encoding a protein folding problem.
        """
        qubit_op = self._qubit_op_builder.build_qubit_op()
        reduced_qubit_op, unused_qubits = remove_unused_qubits(qubit_op)
        self._unused_qubits = unused_qubits
        return reduced_qubit_op

    def interpret(self, raw_result: MinimumEigensolverResult) -> ProteinFoldingResult:
        """
        Interprets the raw algorithm result and returns a ProteinFoldingResult object.

        Args:
            raw_result: A raw result of the protein folding problem.

        Returns:
            A ProteinFoldingResult object that includes the interpreted result.
        """
        best_turn_bitstring = raw_result.best_measurement["bitstring"]
        return ProteinFoldingResult(
            peptide=self._peptide,
            unused_qubits=self._unused_qubits,
            turn_bitstring=best_turn_bitstring,
        )

    @property
    def unused_qubits(self) -> list[int]:
        """Returns the list of indices for qubits in the original problem formulation that were removed during compression."""
        return self._unused_qubits

    @property
    def peptide(self) -> Peptide:
        """Returns the peptide defining the protein subject to the folding problem."""
        return self._peptide
