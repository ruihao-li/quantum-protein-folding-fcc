# Changelog

All notable user-facing changes to this project are documented here.

## Unreleased

### Changed

- Support Qiskit 2.0 through 2.4 using SamplerV2 and EstimatorV2 public APIs.
- Install Qiskit IBM Runtime 0.44 with the core dependencies so the notebook's
  hardware execution mode works after a standalone project installation.
- Replace the deprecated `RealAmplitudes` class with `real_amplitudes()` in the
  workflow, examples, and tests.
- Replace the private `qiskit._accelerate` dependency with NumPy-based Pauli
  term deduplication.
- Accept SamplerV2 primitive results and direct count mappings when interpreting
  folding results.
- Make Ray an optional `ray` installation extra and import it only when the Ray
  parallelizer is selected. The workflow notebook now defaults to an adaptive
  policy: Python multiprocessing on Linux and serial evaluation on macOS or
  Windows, avoiding repeated process-spawn overhead on those platforms.
- Require NumPy 2.x for compatibility with the shared QBioCode and Qiskit
  Machine Learning environment.
- Retain psutil 5.x compatibility with Galaxy's Jupyter resource monitor.
- Add full simulator and IBM Runtime hardware execution modes for both notebook
  workflows, with separate QPU refinement budgets.
- Consolidate all simulator, optimizer, parallelizer, random-seed, and hardware
  settings into a self-contained notebook configuration cell.
- Save simulator-optimized PolyFit and VQEC parameters and dual variables to a
  validated warm-start file for subsequent hardware refinement.
- Allow chance-constrained VQEC to keep a logical ansatz separate from its
  measured, transpiled hardware sampling circuit.

### Compatibility

- Qiskit 2.5 remains intentionally excluded pending a separate compatibility
  validation; this release targets the shared Qiskit 2.2 environment.

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
