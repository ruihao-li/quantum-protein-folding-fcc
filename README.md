# quantum-protein-folding-fcc
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0) [![Supported Python Versions](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/) ![Qiskit](https://img.shields.io/badge/qiskit-2.x-green.svg)


Link to the paper: [Quantum Algorithm for Protein Structure Prediction Using the Face-Centered Cubic Lattice](https://doi.org/10.48550/arXiv.2507.08955)


## Overview

This repository contains the core functionality for performing protein structure prediction with variational quantum algorithms, for a coarse-grained model on the face-centered cubic (FCC) lattice. The implementation is inspired by the existing [quantum protein folding codebase](https://github.com/qiskit-community/quantum-protein-folding) developed by IBM Quantum, based on the article: [Resource-efficient quantum algorithm for protein folding](https://www.nature.com/articles/s41534-021-00368-4). While stylistically similar, our work presents the following key differences:

- The original code is built on the tetrahedral lattice, while this implementation is adapted for the face-centered cubic (FCC) lattice, and can be easily extended to other lattices.
- The `qubit_op_builder` module contains some key differences due to a different formulation of the Hamiltonian for the FCC lattice.
- We include two ways to enforce the non-overlapping constraint without invoking slack variables: the polynomial fit method and a turn-only, chance-constrained VQEC method based on Lagrangian duality.
- In the polynomial fit approach, we use the Sampler primitive. This allows us to accelerate the energy computations with parallelization using [Ray](https://docs.ray.io/en/latest/index.html) and record the best solutions throughout the optimization process.

---

## Turn-only chance-constrained VQEC

The revised VQEC path uses only the compact FCC turn register,
`4 * len(sequence) - 10` qubits. Interaction energies are scored directly from
decoded lattice geometry when `D_mn == 2`; pairwise chance constraints bound
the sampled probability of exact overlap, `D_mn == 0`. Contact ancillas and
optimizer-side postselection are not used.

```python
import numpy as np
from qiskit.circuit.library import real_amplitudes
from qiskit_aer.primitives import SamplerV2 as Sampler

from fcc import build_turn_only_fcc_model
from vqe import ChanceConstrainedVQEC

model = build_turn_only_fcc_model("KLVFFA")
ansatz = real_amplitudes(
    model.num_qubits, reps=1, entanglement="linear"
)

# One explicit overlap-probability limit for each model.constrained_pairs item.
delta_mn = np.full(model.constraint_count, 0.01)
solver = ChanceConstrainedVQEC(
    model,
    ansatz,
    Sampler(default_shots=8192, seed=7),
    constraint_limits=delta_mn,
    shots=8192,
)
result = solver.optimize_primal_dual(
    initial_dual_vars=np.zeros(model.constraint_count),
    max_iter=100,
)
```

The historical `ProteinFoldingProblem.qubit_op()` plus
`olap_constr_ops()` estimator workflow remains available for compatibility,
but it constrains expected distance expressions rather than pairwise overlap
probabilities. See the VQEC section of `workflow_demo.ipynb` for the revised
API.

---

## Quick Start

1. Clone the repository:
   ```bash
   git clone https://github.com/ruihao-li/quantum-protein-folding-fcc.git
   ```
2. Navigate to the project directory and install the required dependencies:
   ```bash
   cd quantum-protein-folding-fcc
   pip install -e .
   ```
   Ray-based parallel processing is optional. Install it with:
   ```bash
   pip install -e ".[ray]"
   ```
3. [Optional] Run the workflow demo notebook to see how it works:
   ```bash
   jupyter notebook workflow_demo.ipynb
   ```

The package supports Qiskit 2.0 through 2.4. Qiskit 2.5 requires NumPy 2,
while this release retains NumPy 1.x compatibility for the broader scientific
Python stack used by the workflow.

---

## Core developers
The quantum components of this codebase are largely developed by [Ruihao Li](https://github.com/ruihao-li) (lir9@ccf.org). The classical exhaustive search methods are developed by [Frank DiFilippo](https://github.com/difilif) (difilif@ccf.org).

Other main contributors:
- [Hakan Doga](https://github.com/hkndoga) (hakandoga@ibm.com)
- [Bryan Raubenolt](https://github.com/thepineapplepirate) (raubenb@ccf.org)
- [Abdullah Ash Saki](https://github.com/ashsaki) (saki@ibm.com)
- [Tomas Radivoyevitch](https://github.com/radivot) (radivot@gmail.com)
- [Daniel Blankenberg](https://github.com/blankenberg) (blanked2@ccf.org)
