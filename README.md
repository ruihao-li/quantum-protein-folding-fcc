# protein-folding-fcc

[![Supported Python Versions](https://img.shields.io/badge/Python-3.10%20%7C%203.10%20%7C%203.12-blue)](https://www.python.org/downloads/) 

## Overview

This repository contains the core functionality for performing protein structure prediction with variational quantum algorithms, for a coarse-grained model on the face-centered cubic (FCC) lattice. The implementation is heavily inspired by the existing [quantum protein folding codebase](https://github.com/qiskit-community/quantum-protein-folding) developed by IBM Quantum, based on the paper: [Resource-efficient quantum algorithm for protein folding](https://www.nature.com/articles/s41534-021-00368-4). While stylistically similar, our work presents the following key differences:

- The original code is built on the tetrahedral lattice, while this implementation is adapted for the face-centered cubic (FCC) lattice, and can be easily extended to other lattices.
- The `qubit_op_builder` module contains some key differences due to a different formulation of the Hamiltonian for the FCC lattice.
- We include two ways to enforce the non-overlapping constraint without invoking slack variables: the polynomial fit and Lagrangian duality methods.
- In running the VQE algorithm with the polynomial fit approach, we use the Sampler primitive. This allows us to accelerate the energy computations with parallelization using `Ray` and record the best solutions throughout the optimization process.

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
3. [Optional] Run the workflow demo notebook to see how it works:
   ```bash
   jupyter notebook workflow_demo.ipynb
   ```

---

## Core developers
The core quantum components of this project are largely developed by [Ruihao Li](https://github.com/ruihao-li) (lir9@ccf.org). The classical exhaustive search methods are developed by [Frank DiFilippo](https://github.com/difilif) (difilif@ccf.org).

Other main contributors:
- [Bryan Raubenolt](https://github.com/thepineapplepirate) (raubenb@ccf.org)
- [Hakan Doga](https://github.com/hkndoga) (hakandoga@ibm.com)
- [Abdullah Ash Saki](https://github.com/ashsaki) (saki@ibm.com)
- [Tomas Radivoyevitch](https://github.com/radivot) (radivot@gmail.com)
