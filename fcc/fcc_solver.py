import time
import numpy as np

from qiskit.primitives import BaseSamplerV2 as BaseSampler
from qiskit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp
from scipy.optimize import minimize
# if you plan on using Ray, uncomment line 9 below, and comment out line 10
# from .measurement_utils import get_cvar_energy, process_counts
from .measurement_utils_multiprocessing import get_cvar_energy, process_counts


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
    ):
        """
        Initialize the solver.

        Args:
            ansatz: The ansatz circuit to use (ansatz circuit must contain
            measurements).
            hamiltonian: The Hamiltonian to use.
            sampler: The sampler to use.
        """
        self.ansatz = ansatz
        self.hamiltonian = hamiltonian
        self.sampler = sampler
        self.cost_trajectory: list[float] = []
        self.global_bitstring_energies: dict[str, float] = {}

    def cost_function(
        self, params: np.ndarray, num_batches: int | None = None, verbose: bool = False
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
        pub_result = primitive_result[0].data
        counts = pub_result.meas.get_counts()  # count occurrences of each bitstring
        toc1 = time.time()
        job_duration = toc1 - tic0

        # 2. Process the counts, calculating the energy of unique bitstrings that have not been processed before
        tic2 = time.time()

        state_wise_energies, prob_energy_pairs = process_counts(
            counts,
            observable=self.hamiltonian,
            num_batches=num_batches,
            global_bitstring_energies=self.global_bitstring_energies,
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
        init_params = (
            np.random.uniform(0, 2 * np.pi, self.ansatz.num_parameters)
            if init_params is None
            else init_params
        )
        optim_result = minimize(
            fun=self.cost_function,
            x0=init_params,
            method=optimizer,
            options={"maxiter": maxiter},
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
            "opt_params": opt_params,
            "cost_trajectory": self.cost_trajectory,
            "top_solutions": top_solutions,
        }
        return final_results
