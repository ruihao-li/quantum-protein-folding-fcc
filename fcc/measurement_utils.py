# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

from itertools import chain
from typing import Iterable, Literal, Mapping
import numpy as np
import psutil
import multiprocessing as mp

try:
    import ray
except ImportError:  # pragma: no cover - exercised only when Ray is unavailable
    ray = None
from qiskit.quantum_info import SparsePauliOp
from qiskit.result import Counts


__all__ = [
    "extract_counts_from_result",
    "normalize_bitstring_mapping",
    "normalize_counts",
    "evaluate_sparsepauli",
    "get_cvar_energy",
    "process_counts",
]

# https://github.com/qiskit-community/qiskit-algorithms/blob/3fb69b3051f39602b1e4c3f0229a849a2ff6b8e4/qiskit_algorithms/minimum_eigensolvers/diagonal_estimator.py#L190-L197
_PARITY = np.array(
    [-1 if bin(i).count("1") % 2 else 1 for i in range(256)], dtype=np.float64
)


def extract_counts_from_result(result_like: object) -> Counts | dict[str, int]:
    """Extract a counts dictionary from a Sampler V2 pub result or its data object.

    Qiskit stores sampled data under the classical-register name. When
    ``measure_all()`` is used the default register is named ``meas``, but user
    circuits may use a different name. This helper first checks the default
    ``meas`` register and then falls back to the first register-like attribute
    exposing ``get_counts()``.
    """
    data = getattr(result_like, "data", result_like)

    if hasattr(data, "get_counts"):
        return data.get_counts()

    default_register = getattr(data, "meas", None)
    if default_register is not None and hasattr(default_register, "get_counts"):
        return default_register.get_counts()

    for attr_name in dir(data):
        if attr_name.startswith("_"):
            continue
        try:
            register_data = getattr(data, attr_name)
        except Exception:
            continue
        if hasattr(register_data, "get_counts"):
            return register_data.get_counts()

    raise TypeError(
        "Could not extract counts from the sampler result. Ensure the circuit contains measurements."
    )


def _normalize_count_bitstring(state: str, num_qubits: int) -> str:
    """Normalize a measurement bitstring to a bare binary string.

    Qiskit count keys may contain register separators such as spaces when a
    circuit uses multiple classical registers. Strip those separators and pad
    leading zeros so the bitstring length matches the observable size.
    """
    normalized = state.replace(" ", "").replace("_", "")
    if set(normalized) - {"0", "1"}:
        raise ValueError(f"Measurement bitstring contains non-binary characters: {state!r}")
    if len(normalized) > num_qubits:
        raise ValueError(
            "Measurement bitstring is longer than the observable qubit count: "
            f"{len(normalized)} > {num_qubits}."
        )
    return normalized.zfill(num_qubits)


def normalize_bitstring_mapping(
    bitstrings: Mapping[str, int | float], num_qubits: int
) -> dict[str, int | float]:
    """Return a mapping keyed by canonical bitstrings while preserving values.

    This is useful for quasi-probability dictionaries, where truncating values to
    integers would destroy the ranking of states.
    """
    normalized_mapping: dict[str, int | float] = {}
    for state, value in bitstrings.items():
        normalized_state = _normalize_count_bitstring(state, num_qubits)
        normalized_mapping[normalized_state] = (
            normalized_mapping.get(normalized_state, 0) + value
        )
    return normalized_mapping


def normalize_counts(
    counts: Counts | dict[str, int], num_qubits: int
) -> dict[str, int]:
    """Return counts keyed by canonical bitstrings.

    Different sampler backends may format count keys differently (for example by
    inserting spaces between classical registers). Canonicalizing the keys in a
    single place keeps downstream caching, reporting, and result interpretation
    consistent.
    """
    normalized_counts: dict[str, int] = {}
    for state, count in counts.items():
        normalized_state = _normalize_count_bitstring(state, num_qubits)
        normalized_counts[normalized_state] = (
            normalized_counts.get(normalized_state, 0) + int(count)
        )
    return normalized_counts


def evaluate_sparsepauli(state: str, observable: SparsePauliOp) -> float:
    """Evaluate a diagonal observable on a measured bitstring.

    Args:
        state (str): Bitstring for which the observable will be evaluated.
        observable (SparsePauliOp): Observable for computing expectation value.

    Returns:
        float: Expectation value.

    Raises:
        ValueError: If ``observable`` contains X or Y terms and is therefore
            not diagonal in the computational basis.
    """
    if np.any(observable.paulis.x):
        raise ValueError(
            "evaluate_sparsepauli only supports diagonal observables composed of I/Z Paulis."
        )

    state = _normalize_count_bitstring(state, observable.num_qubits)
    state_int = int(state, 2)
    packed_uint8 = np.packbits(observable.paulis.z, axis=1, bitorder="little")
    state_bytes = np.frombuffer(
        state_int.to_bytes(packed_uint8.shape[1], "little"), dtype=np.uint8
    )
    reduced = np.bitwise_xor.reduce(packed_uint8 & state_bytes, axis=1)
    coeffs = observable.coeffs
    if coeffs is None:
        raise ValueError("Observable coefficients are not available.")

    return np.real(np.sum(coeffs * _PARITY[reduced]))


def _evaluate_sparsepauli(state: str, observable: SparsePauliOp) -> float:
    """Backward-compatible private alias for :func:`evaluate_sparsepauli`."""
    return evaluate_sparsepauli(state, observable)


# https://github.com/qiskit-community/qiskit-algorithms/blob/3fb69b3051f39602b1e4c3f0229a849a2ff6b8e4/qiskit_algorithms/minimum_eigensolvers/diagonal_estimator.py#L158
# measurements: Iterable[(bitstring probability, expval)]
def get_cvar_energy(
    measurements: Iterable[tuple[float, float]], alpha: float = 0.1
) -> float:
    """Computes Conditional-Value-at-Risk (CVaR) energy.

    Args:
        measurements (Iterable[tuple[float, float]]): An iterable of 2-tuples,
        where tuple contains probability and expectation value of a bitstring
        without explicitly containing the bitstring. The first element of the
        tuple is the probability of the bitstring. The second element of the
        tuple is the expectation value for that bitstring.
        alpha (float): CVaR aggregation.

    Returns:
        float: CVaR energy.

    Raises:
        ValueError: If ``alpha`` is not in the interval ``(0, 1]``.
    """
    if not 0 < alpha <= 1:
        raise ValueError("alpha must be in the interval (0, 1].")

    sorted_measurements = sorted(measurements, key=lambda x: x[1])

    accumulated_percent = 0.0  # once alpha is reached, stop
    cvar = 0.0
    for probability, expval in sorted_measurements:
        cvar += expval * min(probability, alpha - accumulated_percent)
        accumulated_percent += probability
        if accumulated_percent >= alpha:
            break

    return cvar / alpha


def _calculate_batch_energy(
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
    return [
        evaluate_sparsepauli(bitstring, observable)
        for bitstring in conf_bitstring_list
    ]


if ray is not None:
    calculate_batch_energy_ray = ray.remote(_calculate_batch_energy)
else:  # pragma: no cover - exercised only when Ray is unavailable
    calculate_batch_energy_ray = None


def calculate_batch_energy_mp(
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
    return [
        evaluate_sparsepauli(bitstring, observable)
        for bitstring in conf_bitstring_list
    ]


def process_counts(
    counts: Counts | dict[str, int],
    observable: SparsePauliOp,
    num_batches: int | None = None,
    global_bitstring_energies: dict[str, float] | None = None,
    parallelizer: Literal["ray", "python-mp"] = "python-mp",
) -> tuple[dict[str, float], list[tuple[float, float]]]:
    """
    Process a Counts distribution in parallel batches. First, it
    splits all unique states (bitstrings) in the Counts distribution into
    specified number of batches. Then, for each batch of states, it computes the
    energy per unique state. It also converts integer count of a state to
    probability (float) by dividing the count by total number of shots.

    Args:
        counts (Counts): The counts of a quantum circuit run.
        observable (SparsePauliOp): The observable for expectation value
        evaluation.
        num_batches (int, optional): The number of batches to split the
        unique states into. Defaults to None.
        global_bitstring_energies (dict[str, float] | None, optional): The
        global dictionary that keeps track of all bitstrings ever measured.
        Bitstrings that are already in this dictionary will not be
        processed. If None, an empty dictionary is used.
        parallelizer (str, optional): The parallelizer to use. Defaults to
        "python-mp". Options: "ray", "python-mp".

    Returns:
        state_wise_energies (dict[str, float]): Dictionary where keys are unique
        bitstrings (states) from `Counts` and values are the corresponding
        energies of those states.
        prob_energy_pairs (list[tuple[float, float]]): List of 2-tuples. The
        first element of the tuple is the probability of a state. The second
        element of the tuple is the energy of that corresponding state.
    """
    if global_bitstring_energies is None:
        global_bitstring_energies = {}

    if parallelizer not in ("ray", "python-mp"):
        raise ValueError(
            "Unsupported parallelizer. Supported values are 'ray' and 'python-mp'."
        )
    if num_batches is not None and num_batches <= 0:
        raise ValueError("num_batches must be a positive integer or None.")

    counts = normalize_counts(counts, observable.num_qubits)
    total_shots = sum(counts.values())
    if total_shots <= 0:
        raise ValueError("Counts must contain at least one measurement outcome.")
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

    cpu_count = psutil.cpu_count() or 1
    if num_batches is None or num_batches > cpu_count:
        num_batches = cpu_count
    if num_batches > num_unique_states:
        num_batches = num_unique_states
    batch_size = int(np.ceil(num_unique_states / num_batches))

    if parallelizer == "ray":
        if ray is None or calculate_batch_energy_ray is None:
            raise ImportError(
                "parallelizer='ray' requires the optional 'ray' dependency to be installed. "
                "Use parallelizer='python-mp' to avoid the Ray dependency."
            )
        if not ray.is_initialized():
            ray.init(ignore_reinit_error=True, include_dashboard=False)
        observable_id = ray.put(observable)
        doubled_batch_refs = []
        for i in range(0, num_unique_states, batch_size):
            batch_states = unique_states[i : i + batch_size]
            doubled_batch_refs.append(
                calculate_batch_energy_ray.remote(batch_states, observable_id)
            )
        unique_energies = list(chain(*ray.get(doubled_batch_refs)))
        if num_unique_states != len(unique_energies):
            raise RuntimeError(
                "Ray batch processing returned an unexpected number of energy values."
            )
    elif parallelizer == "python-mp":
        with mp.Pool(processes=num_batches) as pool:
            doubled_batch_refs = []
            for i in range(0, num_unique_states, batch_size):
                batched_states = unique_states[i : i + batch_size]
                doubled_batch_refs.append(
                    pool.apply_async(
                        calculate_batch_energy_mp,
                        args=(
                            batched_states,
                            observable,
                        ),
                    )
                )
            unique_energies = list(
                chain.from_iterable(job.get() for job in doubled_batch_refs)
            )
            if num_unique_states != len(unique_energies):
                raise RuntimeError(
                    "Multiprocessing batch processing returned an unexpected number of energy values."
                )

    state_wise_energies: dict[str, float] = {
        state: energy for state, energy in zip(unique_states, unique_energies)
    }
    # Update the global dictionary
    global_bitstring_energies.update(state_wise_energies)
    all_energies = [global_bitstring_energies[state] for state in states]
    prob_energy_pairs: list[tuple[float, float]] = list(zip(state_probs, all_energies))

    return state_wise_energies, prob_energy_pairs
