#include <iostream>
#include <cmath>
#include <format>

double solve_tricomi_tau(double target) {
    double tau = M_PI / 2.0;
    for (int iter = 0; iter < 20; ++iter) {
        double f = tau - std::sin(tau) - target;
        double fp = 1.0 - std::cos(tau);
        double delta = f / fp;
        tau -= delta;
        if (std::abs(delta) < 1e-15) break;
    }
    return tau;
}

double initial_guess_tricomi(int n, int k) {
    int n_half = n / 2;
    double alpha = (n % 2) - 0.5;
    double nu = 4.0 * n_half + 2.0 * alpha + 2.0;

    double target = (4.0 * n_half - 4.0 * k + 3.0) * M_PI / nu;
    std::cout << std::format("  target for tau equation: {}\n", target);

    double tau = solve_tricomi_tau(target);
    std::cout << std::format("  tau = {}\n", tau);

    double sigma = std::cos(tau / 2.0);
    sigma = sigma * sigma;
    std::cout << std::format("  sigma = cos²(tau/2) = {}\n", sigma);

    double x_squared = nu * sigma;
    std::cout << std::format("  x² (before correction) = nu * sigma = {}\n", x_squared);

    if (sigma < 0.9999) {
        double correction = -(1.0 / (3.0 * nu)) * (
            5.0 / (4.0 * (1.0 - sigma) * (1.0 - sigma))
            - 1.0 / (1.0 - sigma)
            - 0.25
        );
        std::cout << std::format("  correction = {}\n", correction);
        x_squared += correction;
    }

    std::cout << std::format("  x² (final) = {}\n", x_squared);
    return std::sqrt(std::max(0.0, x_squared));
}

int main() {
    std::cout << "n=5, k=0 (should give x ≈ 0.9586):\n";
    double x0 = initial_guess_tricomi(5, 0);
    std::cout << std::format("Result: x = {}\n\n", x0);

    std::cout << "n=5, k=1 (should give x ≈ 2.0202):\n";
    double x1 = initial_guess_tricomi(5, 1);
    std::cout << std::format("Result: x = {}\n", x1);

    return 0;
}
