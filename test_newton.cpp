#include <iostream>
#include <cmath>
#include <format>

// Simple test to debug Newton initial guesses
int main() {
    int n = 5;
    int n_half = n / 2;  // = 2
    double alpha = (n % 2) - 0.5;  // = 0.5
    double nu = 4.0 * n_half + 2.0 * alpha + 2.0;  // = 11.0

    std::cout << std::format("n = {}, n_half = {}, alpha = {}, nu = {}\n", 
                            n, n_half, alpha, nu);

    // For k=0, Tricomi should give the smallest positive node
    for (int k = 0; k < n_half; ++k) {
        double target = (4.0 * n_half - 4.0 * k + 3.0) * M_PI / nu;
        std::cout << std::format("k = {}, target = {}\n", k, target);
    }

    return 0;
}
