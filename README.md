# qmutils: Quantum Mechanics Utilities Library

## Overview

qmutils is a C++ library designed to facilitate quantum mechanics calculations
and simulations. It provides a set of tools and abstractions for working with
quantum models, operators, and expressions, making it easier to implement and
study various quantum systems.

## Features

- Quantum operator representations
- Expression building for quantum Hamiltonians
- Common quantum models (e.g., Hubbard, Heisenberg)
- Numerical methods for solving quantum systems
- Utilities for basis state management
- Parallelization support for large-scale computations

## Installation

To install qmutils, follow these steps:

1. Clone the repository:
   ```
   git clone https://github.com/yourusername/qmutils.git
   ```
2. Navigate to the qmutils directory:
   ```
   cd qmutils
   ```
3. Create a build directory and navigate to it:
   ```
   mkdir build && cd build
   ```
4. Run CMake and build the library:
   ```
   cmake ..
   make
   ```

## Quick Start

Here's a simple example to get you started with qmutils:

```cpp
#include <qmutils/model.hpp>
#include <qmutils/expression.hpp>

int main() {
    // Create a Hubbard model
    HubbardChain model(1.0, 0.1, 4.0, 4);  // mu, t, U, chainLength

    // Get the Hamiltonian
    Expression H = model.hamiltonian();

    // Perform calculations...

    return 0;
}
```

## Contributing

We welcome contributions! Please see our [CONTRIBUTING.md](CONTRIBUTING.md) file
for details on how to contribute to qmutils.

## License

qmutils is released under the MIT License. See the [LICENSE](LICENSE) file for
more details.

## Contact

For questions, issues, or suggestions, please open an issue on our [GitHub
repository](https://github.com/keyehzy/qmutils/issues).