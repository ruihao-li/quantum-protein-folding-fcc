# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

"""Builds Pauli operators of a given size."""

import numpy as np
from qiskit.quantum_info import PauliList, SparsePauliOp


def build_full_identity(num_qubits: int) -> SparsePauliOp:
    """
    Builds a full identity operator of a given size.

    Args:
        num_qubits: Number of qubits on which a full identity operator will be
        created.

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
    Assigns predefined values for turn qubits on positions 0, 1, 2, 3, 6, 7 in
    the FCC chain without the loss of generality. Qubits on these positions are
    considered fixed and not subject to optimization. Note that qubits at these
    positions are fixed to 0, which is equivalent to setting Z_i to identity
    operator.

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
    qubit_indices = [0, 1, 2, 3, 6, 7]
    for term in operator:
        table_z = np.copy(term.paulis.z[0])
        table_x = np.copy(term.paulis.x[0])
        for index in qubit_indices:
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
    Finds indices of qubits in a given operator that are equal to an identity
    operator across all terms, i.e., they are irrelevant for the problem.

    Args:
        operator: An operator whose unused qubits shall be removed, e.g., full
        Hamiltonian for the protein folding problem.

    Returns:
        Indices of qubits in the original Hamiltonian that were unused as
        optimization variables.
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
    Calculates reduced Pauli tables by removing qubits that are not used in the
    optimization.

    Args:
        num_qubits: Number of qubits in the original operator.
        table_x: The x array for the symplectic representation of the original
        operator.
        table_z: The z array for the symplectic representation of the original
        operator.
        unused_qubits: List of indices of qubits that are not used in the
        optimization.

    Returns:
        Reduced z and x arrays for the symplectic representation of the
        operator.
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
        unused_qubits: List of indices of qubits that are not used in the
        optimization.

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


def remove_unused_qubits(
    operator: SparsePauliOp, unused_qubit_indices: list[int] = None
) -> tuple[SparsePauliOp, list[int]]:
    """
    Removes qubits in a given operator that are equal to an identity operator
    across all terms, i.e., they are irrelevant for the problem. It makes the
    number of qubits required for encoding the problem smaller or equal.

    Args:
        operator: An operator whose unused qubits shall be removed, e.g., full
        Hamiltonian for the protein folding problem.
        unused_qubit_indices: A list of indices of qubits that are not used in
        the optimization. If None, the function will find them.

    Returns:
        A tuple consisting of the operator compressed to an equivalent one and
        indices of qubits in the original operator that were unused as
        optimization variables.
    """
    if unused_qubit_indices is None:
        unused_qubits = _find_unused_qubits(operator)
    else:
        unused_qubits = unused_qubit_indices
    if isinstance(operator, SparsePauliOp):
        operator_compressed = _compress_sparse_pauli_op(operator, unused_qubits)
        return operator_compressed, unused_qubits
    else:
        return None, None


def compose_IZ_ops(op_1: SparsePauliOp, op_2: SparsePauliOp) -> SparsePauliOp:
    """
    Compose two SparsePauliOps containing only I and Z terms.

    Args:
        op_1: SparsePauliOp containing only I and Z terms.
        op_2: SparsePauliOp containing only I and Z terms.

    Returns:
        SparsePauliOp representing the composed operator.
    """
    assert op_1.num_qubits == op_2.num_qubits, "Mismatched number of qubits"
    num_qubits = op_1.num_qubits

    # Combine Z terms: XOR for Z locations, AND for phase calculation (not needed)
    z1, z2 = op_1.paulis.z, op_2.paulis.z

    # XOR to determine Z locations in the product
    z_combined = np.logical_xor(z1[:, np.newaxis], z2).reshape((-1, num_qubits))
    # Multiply out coefficients
    coeffs = np.multiply.outer(op_1.coeffs, op_2.coeffs).ravel()

    # Simplify the pauli op here
    array = np.packbits(z_combined, axis=1).astype(np.uint16)
    # Find unique Pauli terms
    # Do not depend on Qiskit's private ``qiskit._accelerate`` implementation.
    # NumPy provides the public, stable operation needed here and returns both
    # representative row indices and the inverse map used to combine terms.
    _, indices, inverses = np.unique(
        array, axis=0, return_index=True, return_inverse=True
    )

    if indices.shape[0] == array.shape[0]:
        # No duplicate operator
        pauli_list = PauliList.from_symplectic(z_combined, np.zeros_like(z_combined))
        return SparsePauliOp(pauli_list, coeffs, copy=False)

    coeffs_combined = np.zeros(indices.shape[0], dtype=coeffs.dtype)
    # Combine coefficients of duplicated Pauli terms
    np.add.at(coeffs_combined, inverses, coeffs)
    is_zero = np.isclose(coeffs_combined, 0, atol=1e-8, rtol=1e-5)
    # Check the edge case that we deleted all Paulis
    # In this case we return an identity Pauli with a zero coefficient
    if np.all(is_zero):
        z = np.zeros((1, num_qubits), dtype=bool)
        coeffs = np.array([0j], dtype=coeffs.dtype)
    else:
        non_zero = np.logical_not(is_zero)
        nz_indices = indices[non_zero]
        z = z_combined[nz_indices]
        coeffs = coeffs_combined[non_zero]

    # Create new Pauli list with combined Z locations and phases
    pauli_list = PauliList.from_symplectic(z, np.zeros_like(z))

    return SparsePauliOp(pauli_list, coeffs, ignore_pauli_phase=True, copy=False)
