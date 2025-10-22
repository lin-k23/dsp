#ifndef TRANSFORMS_H
#define TRANSFORMS_H

#include <vector>
#include <complex>

// Type alias
using Complex = std::complex<double>;

// --- Function Declarations ---

/**
 * @brief Computes the DFT using the direct O(N^2) formula.
 */
std::vector<Complex> dft(const std::vector<Complex> &x);

/**
 * @brief Computes the in-place Radix-2 DIT FFT (O(N log N)).
 */
void fft(std::vector<Complex> &a);

#endif // TRANSFORMS_H