
from qiskit.quantum_info import SparsePauliOp, Pauli
from qiskit.quantum_info import commutator
import numpy as np

def mixer_op_x(n, theta: np.pi):
    # Theta = pi for standard X Pauli
    rx_terms = []
    for i in range(n):
        # Create a Pauli string with X on the i-th qubit and I on the others
        labels = ['I'] * n
        labels[i] = 'X'
        pauli_string = ''.join(labels)
        pauli_op = SparsePauliOp(Pauli(pauli_string), coeffs=[theta/2])
        rx_terms.append(pauli_op)
    mixer_hamiltonian = sum(rx_terms, start=SparsePauliOp.from_list([('I'*n, 0)]))
    return mixer_hamiltonian


def truncate_to_two_local(pauli_op):
    # Extract labels and coefficients
    labels = pauli_op.paulis.to_labels()
    coeffs = pauli_op.coeffs
    # Initialize list to collect the filtered terms
    filtered_labels = []
    filtered_coeffs = []
    # Filter terms to only include those with exactly one 'Y' or one 'Y' and one 'Z'
    for label, coeff in zip(labels, coeffs):
        y_count = label.count('Y')
        z_count = label.count('Z')
        # Check for exactly one 'Y', and no 'Z' or exactly one 'Z'
        if (y_count == 1 and z_count == 0) or (y_count == 1 and z_count == 1):
            filtered_labels.append(label)
            filtered_coeffs.append(coeff)
    # Construct the new SparsePauliOp if there are any valid terms
    if filtered_labels:
        truncated_terms = SparsePauliOp.from_list(list(zip(filtered_labels, filtered_coeffs)))
    else:
        # Return an empty operator if no terms meet the criteria
        truncated_terms = SparsePauliOp.from_list([('I' * pauli_op.num_qubits, 0)])
    return truncated_terms


def calculate_commutator(hamiltonian, truncate: bool):
    if truncate:
        # Create the mixer terms first
        mixer_op = mixer_op_x(hamiltonian.num_qubits, np.pi)  # to get regular Pauli X

        # Sets up the version of the adiabatic QC where lambda=0.5 for H_a,
        # partial_H_a is from (1-lambda)H_mixer + lambda * H
        H_a = 0.5 * mixer_op + 0.5 * hamiltonian
        partial_H_a = -mixer_op + hamiltonian

        # Commutator calculated up to second order, l=2
        comm_1 = commutator(H_a, partial_H_a).simplify()  # simplify for cancelling extra terms etc
        comm_2 = commutator(H_a, comm_1).simplify()

        A_2_lambda = comm_1 + comm_2
        A_2_final = A_2_lambda.simplify()

        # Now truncate up to two local terms and extract the ones with Y and YZ
        final_hamiltonian = truncate_to_two_local(A_2_final).simplify()
        return final_hamiltonian

    else:
        # If no truncation is needed, return the simplified version of the second order commutator
        mixer_op = mixer_op_x(hamiltonian.num_qubits, np.pi)  # to get regular Pauli X
        H_a = 0.5 * mixer_op + 0.5 * hamiltonian
        partial_H_a = -mixer_op + hamiltonian

        comm_1 = commutator(H_a, partial_H_a).simplify()
        comm_2 = commutator(H_a, comm_1).simplify()

        A_2_lambda = comm_1 + comm_2
        A_2_final = A_2_lambda.simplify()
        return A_2_final
