#ifndef TRANSFORMS_H
#define TRANSFORMS_H

#include <vector>
#include <complex>

// Define complex number alias
using Complex = std::complex<double>;

/**
 * @brief Compute the Discrete Fourier Transform (DFT) in O(N^2)
 * @param x input time-domain signal
 * @return frequency-domain complex spectrum
 */
std::vector<Complex> dft(const std::vector<Complex> &x);

/**
 * @brief In-place Fast Fourier Transform (FFT) in O(N log N)
 * @param a input/output vector: input is time-domain, output is frequency-domain.
 * @note The vector size N must be a power of two.
 */
void fft(std::vector<Complex> &a);

/**
 * @brief Apply a Hann (Hanning) window to the signal in-place
 * @param x signal vector to modify
 */
void apply_hann_window(std::vector<Complex> &x);

#endif // TRANSFORMS_H
