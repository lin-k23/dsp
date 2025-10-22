#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm> // For std::swap

// Define a type alias for convenience
using Complex = std::complex<double>;

// Define PI
const double PI = std::acos(-1.0);

/**
 * @brief Computes the in-place Radix-2 Decimation-In-Time (DIT) FFT.
 * @param a The input/output vector. Size MUST be a power of 2.
 */
void fft(std::vector<Complex>& a) {
    int N = a.size();
    if (N <= 1) return;

    // --- 0. Get log2(N) ---
    // This determines the number of bits needed for reversal.
    int bits = 0;
    while ((1 << bits) < N) {
        bits++;
    }
    // Basic check that N is a power of 2
    if ((1 << bits) != N) {
        std::cerr << "Error: Input size " << N << " is not a power of 2." << std::endl;
        return;
    }

    // --- 1. Bit-Reversal Permutation ---
    // Reorder the input elements into bit-reversed order.
    for (int i = 0; i < N; ++i) {
        int j = 0;
        // Compute the bit-reversed index 'j' for 'i'
        for (int k = 0; k < bits; ++k) {
            if ((i >> k) & 1) {
                j |= (1 << (bits - 1 - k));
            }
        }
        
        // Swap elements only if i < j to avoid double swaps
        if (i < j) {
            std::swap(a[i], a[j]);
        }
    }

    // --- 2. Butterfly Computations ---
    // Iteratively compute the FFT stages.
    for (int len = 2; len <= N; len <<= 1) { // len = 2, 4, 8, ... N
        int halfLen = len / 2;
        
        // Twiddle factor for this stage: W_len = exp(-j * 2*PI / len)
        Complex W_len_base = std::polar(1.0, -2.0 * PI / len);

        // Iterate over the groups of butterflies
        for (int i = 0; i < N; i += len) {
            Complex W(1.0, 0.0); // W = (W_len_base)^0 = 1

            // Compute the butterflies within this group
            for (int j = 0; j < halfLen; ++j) {
                Complex top = a[i + j];
                Complex bottom = W * a[i + j + halfLen];
                
                a[i + j]         = top + bottom;
                a[i + j + halfLen] = top - bottom;
                
                // Update twiddle factor for next butterfly in group
                W *= W_len_base;
            }
        }
    }
}

// --- Main function for testing ---
int main() {
    // Set N to a power of 2
    int N = 8;
    std::vector<Complex> x(N);

    // Create a sample signal, e.g., x[n] = n
    for (int i = 0; i < N; ++i) {
        x[i] = Complex(i, 0.0);
    }

    std::cout << "Original Signal (x[n]):" << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }

    // --- Compute the FFT ---
    fft(x);
    // --- ---

    std::cout << "\nFFT Result (X[k]):" << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << "X[" << i << "] = " << x[i] << std::endl;
    }

    // Example 2: Impulse signal x = [1, 0, 0, 0]
    // std::vector<Complex> impulse = { {1,0}, {0,0}, {0,0}, {0,0} };
    // fft(impulse);
    // std::cout << "\nFFT of impulse:" << std::endl;
    // for (const auto& val : impulse) std::cout << val << std::endl;
    // (Should print [ (1,0), (1,0), (1,0), (1,0) ])

    return 0;
}