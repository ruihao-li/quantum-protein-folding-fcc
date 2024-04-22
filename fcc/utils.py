"""Builds Pauli operators of a given size."""

from qiskit.quantum_info import SparsePauliOp, PauliList
import numpy as np


def build_full_identity(num_qubits: int) -> SparsePauliOp:
    """
    Builds a full identity operator of a given size.

    Args:
        num_qubits: Number of qubits on which a full identity operator will be created.

    Returns:
        A full identity operator of a given size.
    """
    full_identity = SparsePauliOp("I" * num_qubits)
    return full_identity


def build_pauli_z_op(num_qubits: int, pauli_z_indices: set[int]) -> SparsePauliOp:
    """
    Builds a Pauli operator of a given size with Pauli Z operators on indicated
    positions and identity operators on other positions.

    Args:
        num_qubits: Number of qubits on which a Pauli operator will be created.
        pauli_z_indices: A set of indices in a Pauli operator on which a Pauli Z
        operator shall appear.

    Returns:
        A Pauli operator of a given size with Pauli Z operators on indicated
        positions and identity operators on other positions.
    """
    if 0 in pauli_z_indices:
        operator_str = "Z"
    else:
        operator_str = "I"
    for i in range(1, num_qubits):
        if i in pauli_z_indices:
            operator_str = "Z" + operator_str
        else:
            operator_str = "I" + operator_str
    operator = SparsePauliOp(operator_str)
    return operator


def fix_qubits(operator: SparsePauliOp | int) -> SparsePauliOp | int:
    """
    Assigns predefined values for turn qubits on positions 0, 1, 2, 3, 6, 7 in the FCC chain without the loss of generality. Qubits on these positions are considered fixed and not subject to optimization. Note that qubits at these positions are fixed to 0, which is equivalent to setting Z_i to identity operator.

    Args:
        operator: An operator whose qubits shall be fixed.

    Returns:
        An operator with relevant qubits changed to fixed values.
    """
    if not isinstance(operator, SparsePauliOp):
        return operator
    operator = operator.simplify()
    new_tables_x = []
    new_tables_z = []
    main_bead_indices = [0, 1, 2, 3, 6, 7]
    for term in operator:
        table_z = np.copy(term.paulis.z[0])
        table_x = np.copy(term.paulis.x[0])
        for index in main_bead_indices:
            try:
                # Set Z_i to identity operator means setting table_z[i] to False in the symplectic representation
                table_z[index] = np.bool_(False)
            except IndexError:
                pass
        # table_x is not modified
        new_tables_x.append(table_x)
        new_tables_z.append(table_z)
    new_pauli_table = PauliList.from_symplectic(new_tables_z, new_tables_x)
    operator_updated = SparsePauliOp(data=new_pauli_table, coeffs=operator.coeffs)
    operator_updated = operator_updated.simplify()
    return operator_updated


def _find_unused_qubits(operator: SparsePauliOp) -> list[int]:
    """
    Finds indices of qubits in a given operator that are equal to an identity operator across all terms, i.e., they are irrelevant for the problem.

    Args:
        operator: An operator whose unused qubits shall be removed, e.g., full Hamiltonian for the protein folding problem.

    Returns:
        Indices of qubits in the original Hamiltonian that were unused as optimization variables.
    """
    used_map: dict[int, bool] = {}
    unused_qubits = []
    num_qubits = operator.num_qubits
    # Construct a map of used qubits
    for term in operator:
        table_z = term.paulis.z[0]
        for i in range(num_qubits):
            if table_z[i]:
                used_map[i] = True

    for i in range(num_qubits):
        if i not in used_map:
            unused_qubits.append(i)
    return unused_qubits


def _calc_reduced_pauli_tables(
    num_qubits: int, table_x: np.ndarray, table_z: np.ndarray, unused_qubits: list[int]
) -> tuple[list[bool], list[bool]]:
    """
    Calculates reduced Pauli tables by removing qubits that are not used in the optimization.

    Args:
        num_qubits: Number of qubits in the original operator.
        table_x: The x array for the symplectic representation of the original operator.
        table_z: The z array for the symplectic representation of the original operator.
        unused_qubits: List of indices of qubits that are not used in the optimization.

    Returns:
        Reduced z and x arrays for the symplectic representation of the operator.
    """
    new_table_z = []
    new_table_x = []
    for ind in range(num_qubits):
        if ind not in unused_qubits:
            new_table_z.append(table_z[ind])
            new_table_x.append(table_x[ind])

    return new_table_z, new_table_x


def _compress_sparse_pauli_op(
    operator: SparsePauliOp, unused_qubits: list[int]
) -> SparsePauliOp:
    """
    Compresses a SparsePauliOp by removing unused qubits.

    Args:
        operator: A SparsePauliOp to be compressed.
        unused_qubits: List of indices of qubits that are not used in the optimization.

    Returns:
        A compressed SparsePauliOp.
    """
    new_tables_x = []
    new_tables_z = []
    num_qubits = operator.num_qubits
    for term in operator:
        table_z = term.paulis.z[0]
        table_x = term.paulis.x[0]
        coeffs = term.coeffs[0]
        new_table_z, new_table_x = _calc_reduced_pauli_tables(
            num_qubits, table_x, table_z, unused_qubits
        )
        new_tables_z.append(new_table_z)
        new_tables_x.append(new_table_x)

    new_pauli_table = PauliList.from_symplectic(new_tables_z, new_tables_x)
    operator_compressed = SparsePauliOp(
        data=new_pauli_table, coeffs=operator.coeffs
    ).simplify()
    return operator_compressed


def remove_unused_qubits(operator: SparsePauliOp) -> tuple[SparsePauliOp, list[int]]:
    """
    Removes qubits in a given operator that are equal to an identity operator across all terms, i.e., they are irrelevant for the problem. It makes the number of qubits required for encoding the problem smaller or equal.

    Args:
        operator: An operator whose unused qubits shall be removed, e.g., full Hamiltonian for the protein folding problem.

    Returns:
        A tuple consisting of the operator compressed to an equivalent one and indices of qubits in the original operator that were unused as optimization variables.
    """
    unused_qubits = _find_unused_qubits(operator)
    if isinstance(operator, SparsePauliOp):
        operator_compressed = _compress_sparse_pauli_op(operator, unused_qubits)
        return operator_compressed, unused_qubits
    else:
        return None, None
