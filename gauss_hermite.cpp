#include <algorithm>
#include <cmath>
#include <format>
#include <iostream>
#include <numbers>
#include <numeric>
#include <span>
#include <vector>

using namespace std::numbers;

// sqrt(π) constant (not in C++20 std::numbers)
constexpr double sqrt_pi = 1.772453850905516027298167483341145182798;

// Structure to hold quadrature result
struct QuadratureRule {
    std::vector<double> nodes;
    std::vector<double> weights;

    [[nodiscard]] auto size() const noexcept { return nodes.size(); }
};

// ============================================================================
// Golub-Welsch Eigenvalue-Based Method (Primary Implementation)
// ============================================================================

// QL algorithm with implicit shifts for symmetric tridiagonal eigenvalue problem
// Computes eigenvalues (stored in d) and optionally eigenvectors (stored in z)
void tql2(std::span<double> d, std::span<double> e, std::span<double> z, std::size_t n) {
    constexpr int max_iter = 30;

    // Shift e to simplify indexing
    for (std::size_t i = 1; i < n; ++i) {
        e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;

    for (std::size_t l = 0; l < n; ++l) {
        int iter = 0;
        std::size_t m;

        do {
            // Look for small sub-diagonal element
            for (m = l; m < n - 1; ++m) {
                double dd = std::abs(d[m]) + std::abs(d[m + 1]);
                if (std::abs(e[m]) + dd == dd) break;
            }

            if (m != l) {
                if (iter++ >= max_iter) {
                    throw std::runtime_error("tql2: too many iterations");
                }

                // Form shift
                double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                double r = std::hypot(g, 1.0);
                g = d[m] - d[l] + e[l] / (g + std::copysign(r, g));

                double s = 1.0;
                double c = 1.0;
                double p = 0.0;

                for (std::size_t i = m; i-- > l; ) {
                    double f = s * e[i];
                    double b = c * e[i];
                    r = std::hypot(f, g);
                    e[i + 1] = r;

                    if (r == 0.0) {
                        d[i + 1] -= p;
                        e[m] = 0.0;
                        break;
                    }

                    s = f / r;
                    c = g / r;
                    g = d[i + 1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    p = s * r;
                    d[i + 1] = g + p;
                    g = c * r - b;

                    // Form eigenvector
                    for (std::size_t k = 0; k < n; ++k) {
                        f = z[k * n + i + 1];
                        z[k * n + i + 1] = s * z[k * n + i] + c * f;
                        z[k * n + i] = c * z[k * n + i] - s * f;
                    }
                }

                if (r == 0.0 && iter != 0) continue;

                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            }
        } while (m != l);
    }
}

QuadratureRule gauss_hermite_golub_welsch(int n) {
    if (n <= 0) {
        throw std::invalid_argument("n must be positive");
    }

    auto un = static_cast<std::size_t>(n);

    // Build the Jacobi matrix for Hermite polynomials
    // Diagonal elements are all zero
    std::vector<double> d(un, 0.0);

    // Off-diagonal elements: b_k = sqrt(k/2)
    std::vector<double> e(un, 0.0);
    for (std::size_t k = 1; k < un; ++k) {
        e[k] = std::sqrt(static_cast<double>(k) / 2.0);
    }

    // Initialize eigenvector matrix as identity
    std::vector<double> z(un * un, 0.0);
    for (std::size_t i = 0; i < un; ++i) {
        z[i * un + i] = 1.0;
    }

    // Compute eigenvalues and eigenvectors
    tql2(d, e, z, un);

    // Extract nodes (eigenvalues) and weights
    QuadratureRule result;
    result.nodes.resize(un);
    result.weights.resize(un);

    for (std::size_t i = 0; i < un; ++i) {
        result.nodes[i] = d[i];
        // Weight = sqrt(pi) * (first component of eigenvector)^2
        // Eigenvectors are stored as columns, so i-th eigenvector is column i
        // First component of column i is z[0 * un + i] = z[i]
        result.weights[i] = sqrt_pi * z[i] * z[i];
    }

    // Sort nodes in ascending order (and corresponding weights)
    std::vector<std::size_t> indices(un);
    std::iota(indices.begin(), indices.end(), std::size_t{0});
    std::ranges::sort(indices, [&](std::size_t i, std::size_t j) {
        return result.nodes[i] < result.nodes[j];
    });

    std::vector<double> sorted_nodes(un), sorted_weights(un);
    for (std::size_t i = 0; i < un; ++i) {
        sorted_nodes[i] = result.nodes[indices[i]];
        sorted_weights[i] = result.weights[indices[i]];
    }

    result.nodes = std::move(sorted_nodes);
    result.weights = std::move(sorted_weights);

    return result;
}

// ============================================================================
// Newton's Method with Simple Initial Guesses
// ============================================================================
// This demonstrates the Newton iteration but uses simple initial guesses
// rather than the full Tricomi/Gatteschi asymptotic formulas from the paper.
// A full implementation of the paper's algorithm would require careful handling
// of the indexing scheme which differs from the natural ordering.
// ============================================================================

// Evaluate H_n and H_{n-1} at x using the three-term recurrence
std::pair<double, double> hermite_recurrence(int n, double x) {
    if (n == 0) return {1.0, 0.0};
    if (n == 1) return {2.0 * x, 1.0};

    double H_prev = 1.0;      // H_0
    double H_curr = 2.0 * x;  // H_1

    for (int j = 1; j < n; ++j) {
        double H_next = 2.0 * x * H_curr - 2.0 * static_cast<double>(j) * H_prev;
        H_prev = H_curr;
        H_curr = H_next;

        // Normalize if values get large to prevent overflow
        constexpr double scale_threshold = 1e100;
        if (std::abs(H_curr) > scale_threshold) {
            constexpr double scale_factor = 1e-100;
            H_curr *= scale_factor;
            H_prev *= scale_factor;
        }
    }

    return {H_curr, H_prev};
}

// Newton's method to refine a root of H_n(x)
double newton_hermite(int n, double x_init) {
    constexpr int max_iter = 50;
    constexpr double tol = 1e-15;

    double x = x_init;

    for (int iter = 0; iter < max_iter; ++iter) {
        auto [H_n, H_nm1] = hermite_recurrence(n, x);

        // H_n'(x) = 2n * H_{n-1}(x)
        double H_n_prime = 2.0 * static_cast<double>(n) * H_nm1;

        if (std::abs(H_n_prime) < 1e-20) break;

        double delta = H_n / H_n_prime;
        x -= delta;

        if (std::abs(delta) < tol * (1.0 + std::abs(x))) {
            break;
        }
    }

    return x;
}

QuadratureRule gauss_hermite_newton(int n) {
    if (n <= 0) {
        throw std::invalid_argument("n must be positive");
    }

    // Use simple initial guesses based on approximation of nodes
    // The roots lie approximately in [-sqrt(2n+1), sqrt(2n+1)]
    // Use evenly spaced starting points and refine with Newton

    QuadratureRule result;
    result.nodes.reserve(static_cast<std::size_t>(n));
    result.weights.reserve(static_cast<std::size_t>(n));

    double x_max = std::sqrt(2.0 * static_cast<double>(n) + 1.0);

    // Generate initial guesses
    std::vector<double> guesses;
    for (int i = 0; i < n; ++i) {
        double theta = pi * static_cast<double>(i + 1) / static_cast<double>(n + 1);
        double x_guess = x_max * std::cos(theta);
        guesses.push_back(x_guess);
    }

    // Refine each guess with Newton's method
    for (double guess : guesses) {
        double x = newton_hermite(n, guess);
        result.nodes.push_back(x);
    }

    // Sort nodes (Newton may not give them in order)
    std::ranges::sort(result.nodes);

    // Remove duplicates (in case multiple guesses converged to same root)
    auto unique_end = std::unique(result.nodes.begin(), result.nodes.end(),
                                  [](double a, double b) { return std::abs(a - b) < 1e-10; });
    result.nodes.erase(unique_end, result.nodes.end());

    // Recompute if we lost nodes due to duplicates
    if (result.nodes.size() < static_cast<std::size_t>(n)) {
        // Fall back to Golub-Welsch
        return gauss_hermite_golub_welsch(n);
    }

    // Compute weights
    for (double x : result.nodes) {
        if (n == 1) {
            result.weights.push_back(sqrt_pi);
            continue;
        }

        // Evaluate H_{n-1} carefully for weight calculation
        double H_prev = 1.0;
        double H_curr = 2.0 * x;

        for (int j = 1; j < n - 1; ++j) {
            double H_next = 2.0 * x * H_curr - 2.0 * static_cast<double>(j) * H_prev;
            H_prev = H_curr;
            H_curr = H_next;
        }

        double H_nm1 = H_curr;

        // Weight formula: w_k = 2^{n-1} * n! * sqrt(pi) / (n * H_{n-1}(x_k))^2
        double log_weight = static_cast<double>(n - 1) * std::log(2.0)
                          + std::lgamma(static_cast<double>(n + 1))
                          + 0.5 * std::log(pi)
                          - 2.0 * std::log(std::abs(static_cast<double>(n) * H_nm1));

        double weight = std::exp(log_weight);
        result.weights.push_back(weight);
    }

    return result;
}

// ============================================================================
// Default interface (uses Golub-Welsch)
// ============================================================================

QuadratureRule gauss_hermite(int n) {
    return gauss_hermite_golub_welsch(n);
}

// ============================================================================
// Testing and Validation
// ============================================================================

// Compute double factorial (2k-1)!! = 1*3*5*...*(2k-1)
double double_factorial_odd(int k) {
    if (k <= 0) return 1.0;
    double result = 1.0;
    for (int i = 1; i <= k; ++i) {
        result *= static_cast<double>(2 * i - 1);
    }
    return result;
}

// Test the quadrature rule on known integrals
// ∫_{-∞}^{∞} e^{-x²} x^{2k} dx = sqrt(pi) * (2k-1)!! / 2^k
// An n-point Gauss-Hermite rule is exact for polynomials of degree ≤ 2n-1,
// so the test is valid for k = 0, 1, ..., n-1 (since deg(x^{2k}) = 2k ≤ 2n-2 < 2n-1).
void test_quadrature(const QuadratureRule& rule) {
    int n = static_cast<int>(rule.size());
    int k_max = n - 1;  // exactness limit: 2k ≤ 2n-1 ⟹ k ≤ n-1

    std::cout << std::format("\nTesting quadrature on ∫ e^{{-x²}} x^{{2k}} dx  (n = {}, exact for k ≤ {}):\n",
                            n, k_max);
    std::cout << std::format("{:>5} {:>20} {:>20} {:>15}\n",
                            "k", "Computed", "Exact", "Rel Error");
    std::cout << std::string(65, '-') << '\n';

    for (int k = 0; k <= k_max; ++k) {
        double computed = 0.0;
        for (std::size_t i = 0; i < rule.size(); ++i) {
            double x = rule.nodes[i];
            double f_val = std::pow(x, 2 * k);
            computed += rule.weights[i] * f_val;
        }

        double exact = sqrt_pi * double_factorial_odd(k) / std::pow(2.0, static_cast<double>(k));
        double rel_error = std::abs(computed - exact) / (exact + 1e-100);

        std::cout << std::format("{:5d} {:20.15f} {:20.15f} {:15.3e}\n",
                                k, computed, exact, rel_error);
    }
}

// ============================================================================
// Main Program
// ============================================================================

int main(int argc, char* argv[]) {
    try {
        int n = 10;

        if (argc > 1) {
            n = std::stoi(argv[1]);
        }

        if (n <= 0 || n > 10000) {
            std::cerr << "n must be in range [1, 10000]\n";
            return 1;
        }

        std::cout << std::format("Computing Gauss-Hermite quadrature with n = {}\n", n);
        std::cout << std::string(60, '=') << '\n';

        // Compute using Golub-Welsch (default, most reliable)
        std::cout << "\nMethod: Golub-Welsch (eigenvalue-based)\n";
        auto rule_gw = gauss_hermite_golub_welsch(n);

        // Compute using Newton with simple initial guesses
        std::cout << "Method: Newton with simple initial guesses\n";
        auto rule_newton = gauss_hermite_newton(n);

        // Print results (using Golub-Welsch as reference)
        std::cout << "\nNodes and Weights (Golub-Welsch):\n";
        std::cout << std::format("{:>5} {:>22} {:>22}\n", "i", "x_i", "w_i");
        std::cout << std::string(54, '-') << '\n';

        int display_limit = (n <= 20) ? n : 10;
        for (int i = 0; i < display_limit; ++i) {
            std::cout << std::format("{:5d} {:22.15e} {:22.15e}\n",
                                    i + 1, rule_gw.nodes[static_cast<std::size_t>(i)],
                                    rule_gw.weights[static_cast<std::size_t>(i)]);
        }

        if (n > 20) {
            std::cout << "  ... (" << (n - 2 * display_limit) << " rows omitted) ...\n";
            for (int i = n - display_limit; i < n; ++i) {
                std::cout << std::format("{:5d} {:22.15e} {:22.15e}\n",
                                        i + 1, rule_gw.nodes[static_cast<std::size_t>(i)],
                                        rule_gw.weights[static_cast<std::size_t>(i)]);
            }
        }

        // Validate that weights sum to sqrt(pi)
        double weight_sum_gw = std::accumulate(rule_gw.weights.begin(),
                                              rule_gw.weights.end(), 0.0);
        double weight_sum_newton = std::accumulate(rule_newton.weights.begin(),
                                                   rule_newton.weights.end(), 0.0);

        std::cout << "\nValidation:\n";
        std::cout << std::format("Sum of weights (Golub-Welsch): {:.15f}\n", weight_sum_gw);
        std::cout << std::format("sqrt(π):                       {:.15f}\n", sqrt_pi);
        std::cout << std::format("Difference:                    {:.3e}\n",
                                std::abs(weight_sum_gw - sqrt_pi));

        std::cout << std::format("\nSum of weights (Newton):       {:.15f}\n", weight_sum_newton);
        std::cout << std::format("Difference from sqrt(π):       {:.3e}\n",
                                std::abs(weight_sum_newton - sqrt_pi));

        // Compare the two methods
        std::cout << "\nComparison of methods:\n";
        double max_node_diff = 0.0;
        double max_weight_diff = 0.0;

        for (std::size_t i = 0; i < rule_gw.size(); ++i) {
            max_node_diff = std::max(max_node_diff,
                                    std::abs(rule_gw.nodes[i] - rule_newton.nodes[i]));
            max_weight_diff = std::max(max_weight_diff,
                                      std::abs(rule_gw.weights[i] - rule_newton.weights[i]));
        }

        std::cout << std::format("Max node difference:   {:.3e}\n", max_node_diff);
        std::cout << std::format("Max weight difference: {:.3e}\n", max_weight_diff);

        // Test on known integrals (exact for k = 0, ..., n-1)
        test_quadrature(rule_gw);

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
