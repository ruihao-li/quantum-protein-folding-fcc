# Changelog

All notable user-facing changes to this project are documented here.

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
