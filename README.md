# Gauss-Hermite Quadrature Implementation

This C++20 program computes Gauss-Hermite quadrature nodes and weights for numerical integration of functions with the weight e^{-x²} over the entire real line.

## Overview

Gauss-Hermite quadrature approximates integrals of the form:

```
∫_{-∞}^{∞} e^{-x²} f(x) dx ≈ Σ_{k=1}^{n} w_k f(x_k)
```

where {x_k} are the nodes (roots of the nth Hermite polynomial H_n(x)) and {w_k} are the corresponding weights.

## Implementation

The program implements **two methods**:

### 1. Golub-Welsch Eigenvalue Method (Primary)

This is the gold standard for computing Gauss quadrature rules. It:
- Constructs the symmetric tridiagonal Jacobi matrix for Hermite polynomials
- Computes eigenvalues (nodes) and eigenvectors (for weights) using the QL algorithm
- Provides O(n²) complexity with excellent numerical stability

**Jacobi matrix structure for Hermite polynomials:**
- Diagonal: all zeros
- Off-diagonal: b_k = √(k/2) for k = 1, 2, ..., n-1

### 2. Newton's Method with Initial Guesses

Demonstrates the Newton iteration approach:
- Uses simple trigonometric initial guesses
- Refines roots using Newton's method on H_n(x) = 0
- Evaluates Hermite polynomials via three-term recurrence
- Computes weights using the formula w_k = 2^{n-1} n! √π / [n H_{n-1}(x_k)]²

**Note:** A full implementation of the asymptotic initial guesses from Townsend, Trogdon, and Olver (arXiv:1410.5286v1) would use Tricomi's formula for inner nodes and Gatteschi's formula for outer nodes, but requires careful handling of the paper's indexing scheme.

## Building

```bash
mkdir build
cd build
cmake ..
make
```

## Usage

```bash
./gauss_hermite [n]
```

where `n` is the number of quadrature points (default: 10).

### Examples

```bash
./gauss_hermite 5    # Compute 5-point Gauss-Hermite rule
./gauss_hermite 100  # Compute 100-point rule
```

## Output

The program outputs:
1. **Nodes and weights** computed by both methods
2. **Validation** that weights sum to √π
3. **Comparison** between the two methods
4. **Accuracy tests** on known integrals: ∫ e^{-x²} x^{2k} dx = √π (2k-1)!! / 2^k

### Sample Output (n=10)

```
Nodes and Weights (Golub-Welsch):
    i                    x_i                    w_i
------------------------------------------------------
    1 -3.436159118837739e+00  7.640432855232602e-06
    2 -2.532731674232788e+00  1.343645746781237e-03
    ...
   10  3.436159118837737e+00  7.640432855232649e-06

Validation:
Sum of weights (Golub-Welsch): 1.772453850905515
sqrt(π):                       1.772453850905516
Difference:                    8.882e-16
```

## Features

### C++20 Features Used

- `std::span` for safe array views
- `std::ranges::sort` for expressive sorting
- `std::format` for type-safe formatted output
- `std::numbers::pi` for mathematical constants
- Structured bindings for tuple unpacking
- Concepts for template constraints

### Numerical Techniques

- **QL algorithm** with implicit shifts for symmetric tridiagonal eigenproblems
- **Three-term recurrence** for evaluating Hermite polynomials
- **Overflow prevention** via scaling during polynomial evaluation
- **Log-space computation** for weights to avoid overflow with factorial terms

## Mathematical Background

### Hermite Polynomials

The physicist's Hermite polynomials satisfy:
```
H_{n+1}(x) = 2x H_n(x) - 2n H_{n-1}(x)
H_0(x) = 1
H_1(x) = 2x
```

### Weight Formula

At each node x_k (a root of H_n), the weight is:
```
w_k = (2^{n-1} n! √π) / [n H_{n-1}(x_k)]²
```

### Properties

- Nodes are symmetric: x_k = -x_{n-k+1}
- All nodes lie in (-√(2n+1), √(2n+1))
- Weights sum to √π
- The rule exactly integrates polynomials of degree ≤ 2n-1 with weight e^{-x²}

## References

1. Golub, G. H., & Welsch, J. H. (1969). "Calculation of Gauss quadrature rules." Mathematics of Computation, 23(106), 221-230.

2. Townsend, A., Trogdon, T., & Olver, S. (2014). "Fast computation of Gauss quadrature nodes and weights on the whole real line." arXiv:1410.5286v1.

3. Abramowitz, M., & Stegun, I. A. (1964). "Handbook of Mathematical Functions."

## License

This implementation is provided for educational and research purposes.
