# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

import time
import numpy as np

from qiskit.primitives import BaseSamplerV2 as BaseSampler
from qiskit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp
from scipy.optimize import minimize
from .measurement_utils import (
    extract_counts_from_result,
    get_cvar_energy,
    process_counts,
)


class ProteinSolver:
    """
    Class for solving the protein folding problem using the variational quantum
    eigensolver.
    """

    def __init__(
        self,
        ansatz: QuantumCircuit,
        hamiltonian: SparsePauliOp,
        sampler: BaseSampler,
        parallelizer: str = "python-mp",
    ):
        """
        Initialize the solver.

        Args:
            ansatz: The ansatz circuit to use. If it does not already contain
            measurements, they are added automatically.
            hamiltonian: The Hamiltonian to use.
            sampler: The sampler to use.
            parallelizer: The parallelizer to use. Defaults to "python-mp". Options:
            "ray", "python-mp".
        """
        self.ansatz = ansatz.copy()
        if self.ansatz.num_qubits != hamiltonian.num_qubits:
            raise ValueError(
                "The ansatz and Hamiltonian must act on the same number of qubits: "
                f"{self.ansatz.num_qubits} != {hamiltonian.num_qubits}."
            )
        measured_qubits = {
            self.ansatz.find_bit(instruction.qubits[0]).index
            for instruction in self.ansatz.data
            if instruction.operation.name == "measure"
        }
        if not measured_qubits:
            self.ansatz.measure_all()  # Ensure the ansatz has measurements
        elif measured_qubits != set(range(self.ansatz.num_qubits)):
            raise ValueError(
                "The ansatz must either measure all qubits or contain no measurements. "
                "Partially measured circuits would produce invalid bitstrings for energy evaluation."
            )
        self.hamiltonian = hamiltonian
        self.sampler = sampler
        if parallelizer not in ("ray", "python-mp"):
            raise ValueError(
                "Unsupported parallelizer. Supported values are 'ray' and 'python-mp'."
            )
        self.parallelizer = parallelizer
        self.cost_trajectory: list[float] = []
        self.global_bitstring_energies: dict[str, float] = {}

    def cost_function(
        self,
        params: np.ndarray,
        num_batches: int | None = None,
        verbose: bool = False,
    ) -> float:
        """
        Compute the cost function for the given parameters.

        Args:
            params: The parameters for the ansatz.
            num_batches: The number of batches to split the unique states into.
            Defaults to None. If None, it will be set to the number of CPU
            cores.
            verbose: Whether to print the duration of the job and processing.

        Returns:
            float: The cost.
        """
        # 1. Run the sampler job
        tic0 = time.time()
        pub = (self.ansatz, params)
        job = self.sampler.run(pubs=[pub])
        primitive_result = job.result()
        pub_result = primitive_result[0]
        counts = extract_counts_from_result(pub_result)
        toc1 = time.time()
        job_duration = toc1 - tic0

        # 2. Process the counts, calculating the energy of unique bitstrings that have not been processed before
        tic2 = time.time()
        state_wise_energies, prob_energy_pairs = process_counts(
            counts,
            observable=self.hamiltonian,
            num_batches=num_batches,
            global_bitstring_energies=self.global_bitstring_energies,
            parallelizer=self.parallelizer,
        )

        # 3. Calculate the cost of the measurements
        cost = get_cvar_energy(measurements=prob_energy_pairs, alpha=0.01)
        self.cost_trajectory.append(cost)
        toc2 = time.time()
        process_duration = toc2 - tic2

        if verbose:
            print(f" >> Sampler job took {job_duration:.2f} seconds")
            print(
                f" >> Processing {len(state_wise_energies)} unique bitstrings took {process_duration:.2f} seconds"
            )
        return cost

    def train(
        self,
        optimizer: str,
        maxiter: int,
        num_batches: int | None = None,
        init_params: np.ndarray | None = None,
        num_saved_states: int = 50,
        verbose: bool = False,
    ) -> dict:
        """
        Train the quantum ansatz using the given optimizer.

        Args:
            optimizer: The optimizer to use.
            maxiter: The maximum number of iterations to run the optimizer.
            num_batches: The number of batches to split the unique states into.
            Defaults to None. If None, it will be set to the number of CPU
            cores.
            init_params: The initial parameters for the ansatz. Defaults to None.
            num_saved_states: The number of top states to save. Defaults to 50.
            verbose: Whether to print the duration of the job and processing.

        Returns:
            dict: A dictionary containing the necessary results.
        """
        if not isinstance(optimizer, str) or not optimizer.strip():
            raise ValueError("optimizer must be a non-empty string.")
        if maxiter <= 0:
            raise ValueError("maxiter must be a positive integer.")

        init_params = (
            np.random.uniform(0, 2 * np.pi, self.ansatz.num_parameters)
            if init_params is None
            else np.asarray(init_params, dtype=float)
        )
        if init_params.shape != (self.ansatz.num_parameters,):
            raise ValueError(
                "init_params must be a one-dimensional array whose length matches ansatz.num_parameters."
            )
        if num_saved_states < 0:
            raise ValueError("num_saved_states must be non-negative.")

        optimizer_name = optimizer.upper()
        optimizer_options = {"maxiter": maxiter}
        if optimizer_name == "COBYLA":
            # SciPy internally maps COBYLA's iteration budget to MAXFUN and warns
            # when it is smaller than ``num_parameters + 2``. Normalizing the
            # budget here avoids the warning and makes the effective budget
            # explicit.
            optimizer_options["maxiter"] = max(maxiter, self.ansatz.num_parameters + 2)

        # Each training run should start with a clean optimization history.
        self.cost_trajectory = []
        self.global_bitstring_energies = {}
        if self.ansatz.num_parameters == 0:
            final_cost = self.cost_function(np.array([], dtype=float), num_batches, verbose)
            sorted_solutions = sorted(
                self.global_bitstring_energies.items(), key=lambda x: x[1]
            )
            if num_saved_states > len(sorted_solutions):
                num_saved_states = len(sorted_solutions)
            top_solutions = sorted_solutions[:num_saved_states]
            return {
                "final_cost": final_cost,
                "opt_params": [],
                "init_params": [],
                "cost_trajectory": self.cost_trajectory,
                "top_solutions": top_solutions,
            }
        optim_result = minimize(
            fun=self.cost_function,
            x0=init_params,
            method=optimizer_name,
            options=optimizer_options,
            args=(num_batches, verbose),
        )
        final_cost = optim_result.fun
        opt_params = optim_result.x
        # Sort and save the top N solutions from global_bitstring_energies
        sorted_solutions = sorted(
            self.global_bitstring_energies.items(), key=lambda x: x[1]
        )
        if num_saved_states > len(sorted_solutions):
            num_saved_states = len(sorted_solutions)
        top_solutions = sorted_solutions[:num_saved_states]

        final_results = {
            "final_cost": final_cost,
            "opt_params": opt_params.tolist(),
            "init_params": init_params.tolist(),
            "cost_trajectory": self.cost_trajectory,
            "top_solutions": top_solutions,
        }
        return final_results
