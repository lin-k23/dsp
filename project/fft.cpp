#include "transforms.h"
#include <cmath>
#include <algorithm> // For std::swap
#include <iostream>  // For std::cerr

// Define PI
const double PI = std::acos(-1.0);

void fft(std::vector<Complex>& a) {
    int N = a.size();
    if (N <= 1) return;

    int bits = 0;
    while ((1 << bits) < N) {
        bits++;
    }
    if ((1 << bits) != N) {
        std::cerr << "Error: FFT size " << N << " is not a power of 2." << std::endl;
        return;
    }

    // Bit-Reversal Permutation
    for (int i = 0; i < N; ++i) {
        int j = 0;
        for (int k = 0; k < bits; ++k) {
            if ((i >> k) & 1) {
                j |= (1 << (bits - 1 - k));
            }
        }
        if (i < j) {
            std::swap(a[i], a[j]);
        }
    }

    // Butterfly Computations
    for (int len = 2; len <= N; len <<= 1) {
        int halfLen = len / 2;
        Complex W_len_base = std::polar(1.0, -2.0 * PI / len);
        for (int i = 0; i < N; i += len) {
            Complex W(1.0, 0.0);
            for (int j = 0; j < halfLen; ++j) {
                Complex top = a[i + j];
                Complex bottom = W * a[i + j + halfLen];
                a[i + j]         = top + bottom;
                a[i + j + halfLen] = top - bottom;
                W *= W_len_base;
            }
        }
    }
}