#include "transforms.h"
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip> // for std::setw, std::setprecision, std::fixed
#include <fstream>

// Define PI
const double PI = std::acos(-1.0);

int main()
{
    // --- 1. Analysis parameters ---
    const double Fs = 500.0;       // sampling rate (Hz)
    const int N = (1 << 8);        // number of samples (power of two)
    const double delta_f = Fs / N; // frequency resolution

    std::cout << "Spectrum analysis configuration:" << std::endl;
    std::cout << "Fs = " << Fs << " Hz" << std::endl;
    std::cout << "N = " << N << " samples" << std::endl;
    std::cout << "T_total = " << (N / Fs) << " s" << std::endl;
    std::cout << "Delta_f = " << std::fixed << std::setprecision(5) << delta_f << " Hz" << std::endl;

    // --- 2. Generate time-domain signal ---
    std::vector<Complex> signal(N);
    double t;
    for (int i = 0; i < N; ++i)
    {
        t = (double)i / Fs; // current time (s)
        double s_t = 0.8 * std::sin(2 * PI * 103 * t) +
                     1.0 * std::sin(2 * PI * 107 * t) +
                     0.1 * std::sin(2 * PI * 115 * t);
        signal[i] = Complex(s_t, 0.0);
    }

    // --- 3. Apply Hann window ---
    std::cout << "\nApplying Hann window..." << std::endl;
    apply_hann_window(signal);

    // --- 4. Perform in-place FFT ---
    std::cout << "Performing FFT (N=" << N << ")..." << std::endl;

    // Call in-place FFT; 'signal' will contain the frequency-domain result.
    fft(signal);

    std::cout << "FFT complete." << std::endl;

    // --- 5. Analyze and print results ---
    std::cout << "\n--- Spectrum results (amplitude > 0.01) ---" << std::endl;
    std::cout << "-------------------------------------------------------------" << std::endl;
    std::cout << "Bin (k) | Freq (Hz) | Amplitude (scaled)" << std::endl;
    std::cout << "-------------------------------------------------------------" << std::endl;

    // Only analyze up to N/2 (Nyquist)
    for (int k = 0; k <= N / 2; ++k)
    {
        double freq = k * delta_f;

        // Use 'signal[k]' because fft is in-place
        double magnitude_raw = std::abs(signal[k]);

        // Hann window scaling factor (multiply by 4/N for non-DC bins)
        // DC (k=0) uses factor 2/N
        double amplitude_scaled = (k == 0) ? (magnitude_raw * 2.0 / N) : (magnitude_raw * 4.0 / N);

        // Write output into ../output/spectrum_{N}.csv for matlab plotting
        std::ofstream outfile;
        outfile.open("../output/spectrum_" + std::to_string(N) + "_points.csv", std::ios_base::app);
        outfile << freq << "," << amplitude_scaled << "\n";
        outfile.close();

        // Print only significant peaks
        if (amplitude_scaled > 0.01)
        {
            std::cout << std::setw(7) << k << " | "
                      << std::setw(9) << std::fixed << std::setprecision(3) << freq << " | "
                      << std::setw(16) << std::fixed << std::setprecision(4) << amplitude_scaled
                      << std::endl;
        }
    }

    return 0;
}
