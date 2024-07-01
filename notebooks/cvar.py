from typing import Iterable

from qiskit.quantum_info import SparsePauliOp
from qiskit.result import Counts, sampled_expectation_value

# measurements: Iterable[(bitstring probability, expval)]
def get_cvar_energy(measurements: Iterable[tuple[float, float]], alpha: float = 0.1) -> float:
    # sort by values
    sorted_measurements = sorted(measurements, key=lambda x: x[1])

    accumulated_percent = 0.0  # once alpha is reached, stop
    cvar = 0.0
    for probability, expval in sorted_measurements:
        cvar += expval * min(probability, alpha - accumulated_percent)
        accumulated_percent += probability
        if accumulated_percent >= alpha:
            break

    return cvar / alpha

import time
import numpy as np
import ray
from itertools import chain

_PARITY = np.array([-1 if bin(i).count("1") % 2 else 1 for i in range(256)], dtype=np.float64)

def _evaluate_sparsepauli(state: str, observable: SparsePauliOp) -> float:
    state = int(state, 2)
    packed_uint8 = np.packbits(observable.paulis.z, axis=1, bitorder="little")
    state_bytes = np.frombuffer(state.to_bytes(packed_uint8.shape[1], "little"), dtype=np.uint8)
    reduced = np.bitwise_xor.reduce(packed_uint8 & state_bytes, axis=1)
    return np.real(np.sum(observable.coeffs * _PARITY[reduced]))

# Let's start Ray
ray.init(
    num_cpus=8,
    log_to_driver=False,
    ignore_reinit_error=True,
)

@ray.remote
def calc_expval_parallel(batch_states, oper):
    return [_evaluate_sparsepauli(state, observable=oper) for state in batch_states]

def process_counts(counts, observable: SparsePauliOp):
    total_shots = sum(counts.values())
    state_wise_expval = {}
    prob_expval_list = []
    states = list(counts.keys())
    num_unique_bitstrings = len(states)
    print(f'Number of unique bitstrings {num_unique_bitstrings}')
    state_count = np.array(list(counts.values()))
    state_probs = state_count/total_shots

    num_batches = 8
    batch_size = num_unique_bitstrings // num_batches
    results = []
    for i in range(0, num_unique_bitstrings, batch_size):
        batch_states = states[i:i+batch_size]
        results.append(calc_expval_parallel.remote(batch_states, observable))

    expvals = ray.get(results)
    expvals = list(chain.from_iterable(expvals))
    print(f'Num expvals {len(expvals)}')
    assert num_unique_bitstrings == len(expvals)

    state_wise_expval = {state: expval for state, expval in zip(states, expvals)}
    prob_expval_list = [(prob, expval) for prob, expval in zip(state_probs, expvals)]
    
    return state_wise_expval, prob_expval_list

def process_counts_serial(counts, observable: SparsePauliOp):
    total_shots = sum(counts.values())
    state_wise_expval = {}
    prob_expval_list = []

    for state, count in counts.items():
        prob = count/total_shots
        # much slower
        # expval = sampled_expectation_value(dist={state: 1}, oper=observable)
        
        # faster
        expval = _evaluate_sparsepauli(state, observable=observable)

        state_wise_expval[state] = expval
        prob_expval_list.append((prob, expval))

    return state_wise_expval, prob_expval_list

if __name__ == "__main__":
    pass
