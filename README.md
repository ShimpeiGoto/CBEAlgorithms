# CBEAlgorithms
## Overview
Implementation of Controlled Bond Expansion (CBE) Density Matrix Renormalzation Group (DMRG) and CBE Time-Dependent Variational Principle (TDVP) algorithms based on [ITensor C++ library](https://github.com/ITensor/ITensor). The CBE-DMRG and CBE-TDVP algorithms are introduced in [this PRL paper](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.130.246402) ([arXiv vesion](https://arxiv.org/abs/2207.14712)) and [this paper](https://arxiv.org/abs/2208.10972), respectively. These algorithms are efficient for systems with large local Hilbert space such as bosonic systems.

These codes are written for personal reasons and not well tested. Please use at your own risk.

## Requirement
- ITensor C++ library (developed with version 3.2.0)

## Usage
See dmrg.cc and tdvp.cc for basic usage.
The bahaviors of algorithms are controlled by three paramters
- Cutoff: Used for SVD of matrix product states
- CutoffOrth: Used for SVD of matricies expanded by local MPO
- CutoffExpand : Used for SVD determining expanded linear space

## License
The MIT license