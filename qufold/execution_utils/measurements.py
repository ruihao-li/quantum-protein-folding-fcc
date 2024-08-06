from itertools import chain
from typing import Iterable

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


@ray.remote
def calc_expval(batched_states: list[str], observable: SparsePauliOp) -> list[float]:
    """Evaluates expectation values for many states.

    Args:
        batched_states (list[str]): List of bitstrings (states) to evaluate.
        observable (SparsePauliOp): Observable for expectation value evaluation.

    Returns:
        List of float: Expectation values for bitstrings.
    """
    return [
        _evaluate_sparsepauli(state, observable=observable) for state in batched_states
    ]


def process_counts_parallel(
    counts: Counts | dict[str, int],
    observable: SparsePauliOp,
    num_batches: int | None = None,
) -> tuple[dict[str, float], list[tuple[float, float]]]:
    # TODO: split into multiple functions as it is doing many things?
    """Process a Counts distribution in parallel batches.

    First, it splits all unique states (bitstrings) in the Counts distribution into
    specifed number of batches. Then, for each batch of states, it computes expectation
    value per unique state. It also converts integer count of a state to probability
    (float) by dividing the count by total number of shots.

    Finally, it organizes the processed data into two seperate data structure for ease.

    Args:
        counts (Counts | dict[str, int]): Distribution from a Sampler job.
        observable (SparsePauliOp): Observable for expectation value evaluation.
        num_batches (int): Number of batches to be processed in parallel by `ray`.

    Returns:
        state_wise_expval (dict[str, float]): Dictionary where keys are unique
            bitstrings (states) from `Counts` and values are corresponding expectation
            value of the respective state.

        prob_expval_pairs (list[tuple[float, float]]): List of 2-tuples.
            The first element of the tuple is the probability of a state.
            The seconds element of the tuple is the expecation value of that
            corresponding state.
    """
    total_shots = sum(counts.values())
    states = list(counts.keys())
    num_unique_states = len(counts)

    state_count = np.array(list(counts.values()))
    state_probs = state_count / total_shots

    if num_batches is None:
        num_batches = psutil.cpu_count()

    batch_size = num_unique_states // num_batches
    results = []
    for i in range(0, num_unique_states, batch_size):
        batched_states = states[i : i + batch_size]
        results.append(calc_expval.remote(batched_states, observable))

    expvals = list(chain.from_iterable(ray.get(results)))

    assert num_unique_states == len(expvals)

    state_wise_expval: dict[str, float] = {
        state: expval for state, expval in zip(states, expvals) if expval < 50
    }
    prob_expval_pairs: list[tuple[float, float]] = list(zip(state_probs, expvals))

    return state_wise_expval, prob_expval_pairs


def process_counts_serial(
    counts: Counts | dict[str, int], observable: SparsePauliOp
) -> tuple[dict[str, float], list[tuple[float, float]]]:
    # TODO: split into multiple functions as it is doing many things?
    """Process bitstrings (states) from a Counts distribution serailly.

    For each state (bitstring) in Counts, it computes expectation
    value. It also converts integer count of a state to probability
    (float) by dividing the count by total number of shots.

    Finally, it organizes the processed data into two seperate data structure for ease.

    Args:
        counts (Counts | dict[str, int]): Distribution from a Sampler job.
        observable (SparsePauliOp): Observable for expectation value evaluation.
        num_batches (int): Number of batches to be processed in parallel by `ray`.

    Returns:
        state_wise_expval (dict[str, float]): Dictionary where keys are unique
            bitstrings (states) from `Counts` and values are corresponding expectation
            value of the respective state.

        prob_expval_pairs (list[tuple[float, float]]): List of 2-tuples.
            The first element of the tuple is the probability of a state.
            The seconds element of the tuple is the expecation value of that
            corresponding state.
    """
    total_shots = sum(counts.values())
    state_wise_expval = {}
    prob_expval_list = []

    for state, count in counts.items():
        prob = count / total_shots
        expval = _evaluate_sparsepauli(state, observable=observable)

        state_wise_expval[state] = expval
        prob_expval_list.append((prob, expval))

    return state_wise_expval, prob_expval_list
