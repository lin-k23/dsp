#include <iostream>
#include <vector>
#include <complex>
#include <chrono>  // For timing
#include <iomanip> // For std::setprecision
#include <cmath>   // For sin

// Include the header file with our function declarations
#include "transforms.h"

// Main function for benchmarking
int main()
{
    // --- Configuration --
    int N = 1 << 8;                    // 4096 points
    double epsilon = 1e-6;             // Tolerance for verifications
    const double PI = std::acos(-1.0); // main() needs PI for the test signal

    std::cout << "Benchmarking DFT vs. FFT for N = " << N << " points." << std::endl;

    // --- Create Input Signals ---
    std::vector<Complex> signal_dft(N);
    std::vector<Complex> signal_fft(N);
    for (int i = 0; i < N; ++i)
    {
        signal_dft[i] = Complex(std::sin(2.0 * PI * i / 10.0), 0.0); // Sample sine wave
        signal_fft[i] = signal_dft[i];
    }

    // --- Benchmark Direct DFT (O(N^2)) ---
    std::cout << "\nStarting Direct DFT (O(N^2))... " << std::flush;
    auto start_dft = std::chrono::high_resolution_clock::now();

    std::vector<Complex> result_dft = dft(signal_dft);

    auto end_dft = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> time_dft = end_dft - start_dft;
    std::cout << "Done." << std::endl;

    // --- Benchmark FFT (O(N log N)) ---
    std::cout << "Starting FFT (O(N log N))...   " << std::flush;
    auto start_fft = std::chrono::high_resolution_clock::now();

    fft(signal_fft); // In-place

    auto end_fft = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> time_fft = end_fft - start_fft;
    std::cout << "Done." << std::endl;

    // --- Print Timings ---
    std::cout << "\n--- Results ---" << std::endl;
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Direct DFT (O(N^2)) took: " << time_dft.count() << " ms" << std::endl;
    std::cout << "Radix-2 FFT (O(N log N)) took: " << time_fft.count() << " ms" << std::endl;
    std::cout << "-----------------" << std::endl;
    std::cout << "FFT was " << (time_dft.count() / time_fft.count()) << " times faster." << std::endl;

    // --- Verify Correctness ---
    bool all_match = true;
    for (int i = 0; i < N; ++i)
    {
        if (std::abs(result_dft[i] - signal_fft[i]) > epsilon)
        {
            all_match = false;
            std::cout << "Mismatch at index " << i << "!" << std::endl;
            std::cout << "  DFT: " << result_dft[i] << std::endl;
            std::cout << "  FFT: " << signal_fft[i] << std::endl;
            break;
        }
    }

    if (all_match)
    {
        std::cout << "\nVerification: SUCCESS. Outputs match." << std::endl;
    }
    else
    {
        std::cout << "\nVerification: FAILURE. Outputs do not match." << std::endl;
    }

    return 0;
}