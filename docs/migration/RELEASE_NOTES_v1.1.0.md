# quantum-protein-folding-fcc 1.1.0

Version 1.1.0 adds a compact, chance-constrained VQEC workflow for FCC protein structure prediction.

The revised model uses only `4N - 10` turn qubits. It evaluates interaction energy directly from sampled lattice geometries, rewards contacts only when `D_mn == 2`, and treats exact overlaps (`D_mn == 0`) as pairwise probability constraints with caller-supplied limits. No interaction ancillas or optimizer-side postselection are used.

The new `fcc.build_turn_only_fcc_model` and `vqe.ChanceConstrainedVQEC` APIs reuse the package's FCC decoder, MJ interaction model, turn penalties, and existing perturbed primal-dual equations. The demo
notebook and migration notes show the complete workflow.

This is an additive release. Existing public imports, the historical symbolic VQEC estimator path, and PolyFit behavior remain available. The package version is aligned with the existing `v1.0.0` release line; it advances from that public release to `1.1.0` despite the stale `0.1.0` metadata in the previous source tree.

The source-distribution manifest also includes the requirements file consumed by `setup.py`, so the standard PEP 517 wheel-from-sdist build now completes.
