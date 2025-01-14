from itertools import chain
from typing import Iterable

import math
import multiprocessing as mp
import numpy as np
import psutil

import ray
from qiskit.quantum_info import SparsePauliOp
from qiskit.result import Counts

# https://github.com/qiskit-community/qiskit-algorithms/blob/3fb69b3051f39602b1e4c3f0229a849a2ff6b8e4/qiskit_algorithms/minimum_eigensolvers/diagonal_estimator.py#L190-L197
_PARITY = np.array(
    [-1 if bin(i).count("1") % 2 else 1 for i in range(256)], dtype=np.float64
)


def _evaluate_sparsepauli(state: str, observable: SparsePauliOp) -> float:
    """Evaluates an observable for a given bitstring.

    Args:
        state (str): Bitstring for which the observable will be evaluated.
        observable (SparsePauliOp): Observable for computing expectation value.

    Returns:
        float: Expectation value.
    """
    state = int(state, 2)
    packed_uint8 = np.packbits(observable.paulis.z, axis=1, bitorder="little")
    state_bytes = np.frombuffer(
        state.to_bytes(packed_uint8.shape[1], "little"), dtype=np.uint8
    )
    reduced = np.bitwise_xor.reduce(packed_uint8 & state_bytes, axis=1)

    return np.real(np.sum(observable.coeffs * _PARITY[reduced]))


# https://github.com/qiskit-community/qiskit-algorithms/blob/3fb69b3051f39602b1e4c3f0229a849a2ff6b8e4/qiskit_algorithms/minimum_eigensolvers/diagonal_estimator.py#L158
# measurements: Iterable[(bitstring probability, expval)]
def get_cvar_energy(
    measurements: Iterable[tuple[float, float]], alpha: float = 0.1
) -> float:
    """Computes Conditional-Value-at-Risk (CVaR) energy.

    Args:
        measurements (Iterable[tuple[float, float]]): An iterable of 2-tuples, where
            tuple contains probability and expectation value of a bitstring without
            explicitly containing the bitstring. The first element of the tuple is the
            probability of the bitstring. The second element of the tuple is the
            expectation value for that bitstring.
        alpha (float): CVaR aggregation.

    Returns:
        float: CVaR energy.
    """
    sorted_measurements = sorted(measurements, key=lambda x: x[1])

    accumulated_percent = 0.0  # once alpha is reached, stop
    cvar = 0.0
    for probability, expval in sorted_measurements:
        cvar += expval * min(probability, alpha - accumulated_percent)
        accumulated_percent += probability
        if accumulated_percent >= alpha:
            break

    return cvar / alpha


#@ray.remote
def calculate_batch_energy(
    conf_bitstring_list: list[str],
    observable: SparsePauliOp,
) -> list[float]:
    """
    Calculate the energy of a batch of bitstrings.

    Args:
        conf_bitstring_list (list[str]): List of conformation bitstrings.
        observable (SparsePauliOp): The observable to evaluate the energy.

    Returns:
        list[float]: List of energies of the batch of bitstrings.
    """
    # result = []
    # with mp.Pool() as pool:
    #     # first argument below must be the function you're trying to parallelize
    #     # the second argument must always be a list that you're trying to iterate over
    #     new_conf_bitstring_list = [(bitstring, observable) for bitstring in conf_bitstring_list]
    #     calc_batch_energy_result = pool.starmap(_evaluate_sparsepauli, new_conf_bitstring_list) # starmap ensures correct order is preserved
    # return calc_batch_energy_result
    
    return [
        _evaluate_sparsepauli(bitstring, observable)
        for bitstring in conf_bitstring_list
    ]


def process_counts(
    counts: Counts | dict[str, int],
    observable: SparsePauliOp,
    num_batches: int | None = None,
    global_bitstring_energies: dict[str, float] = {},
) -> tuple[dict[str, float], list[tuple[float, float]]]:
    """
    Process a Counts distribution in parallel batches. First, it splits all
    unique states (bitstrings) in the Counts distribution into specified number
    of batches. Then, for each batch of states, it computes the energy per
    unique state. It also converts integer count of a state to probability
    (float) by dividing the count by total number of shots.

    Args:
        counts (Counts): The counts of a quantum circuit run.
        observable (SparsePauliOp): The observable for expectation value evaluation.
        num_batches (int, optional): The number of batches to split the unique
        states into. Defaults to None.
        global_bitstring_energies (dict[str, float], optional): The global
        dictionary that keeps track of all bitstrings ever measured. Bitstrings
        that are already in this dictionary will not be processed. Defaults to
        empty dictionary.

    Returns:
        state_wise_energies (dict[str, float]): Dictionary where keys are unique
        bitstrings (states) from `Counts` (with energy below a certain value)
        and values are corresponding energy of the respective state.
        prob_energy_pairs (list[tuple[float, float]]): List of 2-tuples. The
        first element of the tuple is the probability of a state. The seconds
        element of the tuple is the energy of that corresponding state.
    """
    total_shots = sum(counts.values())
    states = list(counts.keys())
    state_counts = np.array(list(counts.values()))
    state_probs = state_counts / total_shots
    # Get the unique states that are not already in the global dictionary
    unique_states = [
        state for state in states if state not in global_bitstring_energies
    ]
    num_unique_states = len(unique_states)

    if num_unique_states == 0:
        all_energies = [global_bitstring_energies[state] for state in states]
        prob_energy_pairs: list[tuple[float, float]] = list(
            zip(state_probs, all_energies)
        )
        return {}, prob_energy_pairs

    if num_batches is None or num_batches > psutil.cpu_count():
        num_batches = psutil.cpu_count()
    if num_batches > num_unique_states:
        num_batches = num_unique_states
    batch_size = num_unique_states // num_batches
    
    # start parallelizing
    with mp.Pool(processes=num_batches) as pool:
        doubled_batch_refs = [] 
        for i in range(0, num_unique_states, batch_size):
            batched_states = states[i : i + batch_size]
            doubled_batch_refs.append(pool.apply_async(calculate_batch_energy, args=(batched_states, observable,)))

            unique_energies = list(chain.from_iterable([job.get() for job in doubled_batch_refs]))

    assert num_unique_states == len(unique_energies)

    state_wise_energies: dict[str, float] = {
        state: energy for state, energy in zip(unique_states, unique_energies)
    }
    # Update the global dictionary
    global_bitstring_energies.update(state_wise_energies)
    all_energies = [global_bitstring_energies[state] for state in states]
    prob_energy_pairs: list[tuple[float, float]] = list(zip(state_probs, all_energies))

    return state_wise_energies, prob_energy_pairs
