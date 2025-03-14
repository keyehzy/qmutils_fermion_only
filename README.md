# qmutils

[![Passing](https://github.com/keyehzy/qmutils/workflows/CI/badge.svg)]()
[![License](https://img.shields.io/badge/license-MIT-blue)]()

qmutils is a fast, modern C++ library for quantum mechanics calculations. It provides an intuitive interface for quantum operator algebra, including efficient normal ordering and matrix element computation.

## Quick Start

### Installation

```bash
git clone https://github.com/organization/qmutils.git
cd qmutils
mkdir build && cd build
cmake ..
make
```

### Basic Usage

```cpp
#include "qmutils/basis.h"
#include "qmutils/expression.h"
#include "qmutils/matrix_elements.h"
#include "qmutils/operator.h"

using namespace qmutils;

// Create a simple Hamiltonian
const float t = 1.0f;
const float u = 0.5f;
Expression H;
H += t * Term::hopping(0, 1, Operator::Spin::Up);   // Hopping term
H += u * Term::density(Operator::Spin::Up, 0);      // Density term

// Work with basis states
Basis basis(4, 2, /*Sz=*/0);  // 4 orbitals, 2 particles
auto matrix = compute_matrix_elements<SpMat_cf>(basis, H);
```

See [examples](examples/) for more complex usages.

## Requirements

- C++17 compiler (GCC ≥ 7.0, Clang ≥ 5.0, or MSVC ≥ 19.14)
- CMake ≥ 3.10

Optional:
- Google Test (for testing)

## Features
- Creation/annihilation operators with spin and orbital indices.
- Only Fermionic support.
- Efficient symbolic normal ordering with caching.
- Basis state management.
- Matrix element computation (dense and sparse matrix support).
- Utilities for multi-dimensional lattice indexing.


## Documentation

- [Examples](examples/)
- [Docs](docs/)

## License

MIT License - see [LICENSE](LICENSE) for details.

## Citing qmutils

If you use qmutils in your research, please cite:

```bibtex
@software{qmutils2024,
  author = {de Sousa, M.},
  title = {qmutils: A C++ Library for Quantum Many-Body Calculations},
  year = {2024},
  url = {https://github.com/keyehzy/qmutils}
}
```
