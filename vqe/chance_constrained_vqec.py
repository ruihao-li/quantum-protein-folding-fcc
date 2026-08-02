# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

"""Sample-scored pairwise chance constraints for the existing PDP optimizer."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Protocol

import numpy as np
from qiskit import QuantumCircuit
from qiskit.primitives import BaseSamplerV2

from .vqec_optimization import PerturbedPrimalDualOpt


PARAMETER_SHIFT = np.pi / 2.0
DEFAULT_SHOTS = 8192


class SampledChanceConstrainedProblem(Protocol):
    """Interface required by :class:`ChanceConstrainedVQEC`."""

    num_qubits: int
    objective_scale: float

    @property
    def constraint_count(self) -> int:
        """Return the number of pairwise violation primitives."""

    def primitive_values(
        self,
        counts: Mapping[str, int | float]
        | Sequence[Mapping[str, int | float]],
    ) -> np.ndarray:
        """Return objective and violation expectations for each counts mapping."""


class ChanceConstrainedVQEC(PerturbedPrimalDualOpt):
    """Run the existing PDP equations on sample-scored diagonal primitives.

    The ansatz must be an unmeasured circuit whose parameters each obey the
    standard ``+/- pi/2`` shift rule, such as a decomposed ``RealAmplitudes``
    circuit. Every sampled bitstring is passed to the problem without filtering
    or postselection.
    """

    def __init__(
        self,
        problem: SampledChanceConstrainedProblem,
        ansatz: QuantumCircuit,
        sampler: BaseSamplerV2,
        *,
        constraint_limits: Sequence[float] | np.ndarray,
        shots: int = DEFAULT_SHOTS,
    ) -> None:
        if not isinstance(ansatz, QuantumCircuit):
            raise TypeError("ansatz must be a QuantumCircuit")
        if ansatz.num_qubits != int(problem.num_qubits):
            raise ValueError(
                "Ansatz and chance-constrained problem use different qubit counts"
            )
        if ansatz.num_clbits:
            raise ValueError("ChanceConstrainedVQEC expects an unmeasured ansatz")
        if ansatz.num_parameters < 1:
            raise ValueError("The chance-constrained ansatz must be parameterized")
        if not hasattr(sampler, "run"):
            raise TypeError("sampler must implement the SamplerV2 run interface")

        limits = np.asarray(constraint_limits, dtype=float).copy()
        expected_shape = (int(problem.constraint_count),)
        if limits.shape != expected_shape:
            raise ValueError(
                f"constraint_limits has shape {limits.shape}; expected {expected_shape}"
            )
        if not np.isfinite(limits).all() or np.any((limits < 0.0) | (limits > 1.0)):
            raise ValueError("constraint_limits must be finite and lie in [0, 1]")
        limits.setflags(write=False)

        shot_count = int(shots)
        if shot_count < 1 or shot_count != shots:
            raise ValueError("shots must be a positive integer")
        if not math.isfinite(float(problem.objective_scale)) or float(
            problem.objective_scale
        ) <= 0.0:
            raise ValueError("problem.objective_scale must be finite and positive")

        self.problem = problem
        self.constraint_limits = limits
        self.shots = shot_count
        self._sampler = sampler
        self._ansatz = ansatz
        self._measured_ansatz = ansatz.copy()
        self._measured_ansatz.measure_all()

        # The parent update loop uses only the ansatz parameter count and the
        # constraint collection length once these two hooks are overridden.
        self._qubit_op = None
        self._constr_ops = [None] * int(problem.constraint_count)
        self._estimator = None
        self._gradient = None

    def shifted_parameters(
        self, parameters: np.ndarray, *, shift: float = PARAMETER_SHIFT
    ) -> np.ndarray:
        """Return batched plus/minus shifts in circuit-parameter order."""

        values = self._parameter_vector(parameters)
        shift_value = float(shift)
        if not math.isfinite(shift_value):
            raise ValueError("shift must be finite")
        plus = np.tile(values, (self._ansatz.num_parameters, 1))
        minus = plus.copy()
        diagonal = np.arange(self._ansatz.num_parameters)
        plus[diagonal, diagonal] += shift_value
        minus[diagonal, diagonal] -= shift_value
        return np.vstack((plus, minus))

    def sample_counts(
        self, parameters: np.ndarray
    ) -> tuple[Mapping[str, int], ...]:
        """Sample one parameter vector or a batch using SamplerV2."""

        batch = np.atleast_2d(np.asarray(parameters, dtype=float))
        if batch.shape[1:] != (self._ansatz.num_parameters,):
            raise ValueError(
                f"Parameter batch has shape {batch.shape}; expected (*, {self._ansatz.num_parameters})"
            )
        if not np.isfinite(batch).all():
            raise ValueError("Parameter batch contains non-finite values")

        pubs = [(self._measured_ansatz, row) for row in batch]
        result = self._sampler.run(pubs, shots=self.shots).result()
        sampled = []
        for pub_result in result:
            data = pub_result.data
            try:
                counts = data.meas.get_counts()
            except AttributeError as exc:
                raise ValueError(
                    "Sampler result does not contain the expected 'meas' register"
                ) from exc
            normalized = {
                str(bitstring).replace(" ", ""): int(count)
                for bitstring, count in counts.items()
                if int(count) > 0
            }
            if sum(normalized.values()) != self.shots:
                raise ValueError("Sampler counts do not sum to the configured shots")
            sampled.append(normalized)
        if len(sampled) != len(batch):
            raise ValueError("Sampler returned an unexpected number of results")
        return tuple(sampled)

    def primitive_values(self, parameters: np.ndarray) -> np.ndarray:
        """Sample and return normalized objective plus overlap probabilities."""

        return self._primitive_rows(self.sample_counts(parameters))

    def primitive_gradients(self, parameters: np.ndarray) -> np.ndarray:
        """Estimate all primitive gradients with the parameter-shift rule."""

        primitives = self.primitive_values(self.shifted_parameters(parameters))
        count = self._ansatz.num_parameters
        return 0.5 * (primitives[:count] - primitives[count:]).T

    def _evaluate_primitives(self, params: np.ndarray) -> tuple[float, np.ndarray]:
        primitive = self.primitive_values(self._parameter_vector(params))[0]
        residuals = primitive[1:] - self.constraint_limits
        return float(primitive[0]), residuals

    def _evaluate_gradients(
        self, params: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        gradients = self.primitive_gradients(self._parameter_vector(params))
        objective_gradient = gradients[0][None, :]
        constraint_gradients = gradients[1:, None, :]
        return objective_gradient, constraint_gradients

    def _parameter_vector(self, parameters: np.ndarray) -> np.ndarray:
        values = np.asarray(parameters, dtype=float)
        if values.shape == (1, self._ansatz.num_parameters):
            values = values[0]
        if values.shape != (self._ansatz.num_parameters,):
            raise ValueError(
                f"Parameters have shape {values.shape}; expected ({self._ansatz.num_parameters},)"
            )
        if not np.isfinite(values).all():
            raise ValueError("Parameters contain non-finite values")
        return values

    def _primitive_rows(
        self, counts: Sequence[Mapping[str, int | float]]
    ) -> np.ndarray:
        values = np.asarray(self.problem.primitive_values(counts), dtype=float)
        expected = (len(counts), 1 + int(self.problem.constraint_count))
        if values.shape != expected or not np.isfinite(values).all():
            raise ValueError(
                f"Problem primitive_values returned {values.shape}; expected {expected}"
            )
        return values
