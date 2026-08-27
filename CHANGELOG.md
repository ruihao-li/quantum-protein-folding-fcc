# Changelog

All notable user-facing changes to this project are documented here.

## Unreleased

### Changed

- Support Qiskit 2.0 through 2.4 using SamplerV2 and EstimatorV2 public APIs.
- Replace the deprecated `RealAmplitudes` class with `real_amplitudes()` in the
  workflow, examples, and tests.
- Replace the private `qiskit._accelerate` dependency with NumPy-based Pauli
  term deduplication.
- Accept SamplerV2 primitive results and direct count mappings when interpreting
  folding results.
- Make Ray an optional `ray` installation extra and import it only when the Ray
  parallelizer is selected. The workflow notebook now defaults to portable
  Python multiprocessing.
- Relax exact dependency pins while retaining NumPy 1.x compatibility for the
  broader scientific Python stack.
- Add an opt-in, low-shot IBM Runtime SamplerV2 hardware smoke test to the
  workflow notebook.
- Consolidate all simulator, optimizer, parallelizer, random-seed, and hardware
  settings into a self-contained notebook configuration cell.

### Compatibility

- Qiskit 2.5 is intentionally excluded because it requires NumPy 2; this
  release retains compatibility with the shared NumPy 1.x environment.

## 1.1.0

### Added

- A compact FCC turn-only model using `4N - 10` qubits and direct sample-wise interaction scoring at squared lattice distance `D_mn == 2`.
- Pairwise chance constraints for exact overlaps, with an explicit `constraint_limits` value for each generated residue pair.
- A SamplerV2-backed `ChanceConstrainedVQEC` adapter that reuses the existing perturbed primal-dual update equations and parameter-shift workflow.
- Shared compact-turn decoding, coordinate, encoding, and pair-enumeration helpers.
- Focused geometry, primitive, gradient, optimizer-contract, and PolyFit regression tests.

### Changed

- The VQEC section of `workflow_demo.ipynb` now demonstrates the revised turn-only chance-constrained formulation.
- Package metadata now follows the existing public tag line at version 1.1.0, correcting the stale `0.1.0` value left in `setup.py` after release 1.0.0.

### Compatibility

- Existing public imports and the historical `ProteinFoldingProblem.qubit_op()` / `olap_constr_ops()` estimator workflow remain available.
- The legacy full-Hamiltonian and PolyFit path keeps eager geometry-map and contact-ancilla construction by default.
