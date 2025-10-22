#include "transforms.h"
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip> // for std::setw, std::setprecision, std::fixed
#include <fstream>
#include <string> // for std::to_string

// Define PI
const double PI = std::acos(-1.0);

// Helper function to create the output filename
std::string get_filename(const std::string &prefix, int N)
{
    // Note: The user's original file path was "../output/..."
    // We will use "./data/" as requested in a previous prompt.
    return "../output/" + prefix + "_" + std::to_string(N) + "_points.csv";
}

int main()
{
    // --- 1. Analysis parameters ---
    const double Fs = 500.0;       // sampling rate (Hz)
    const int N = (1 << 9);        // number of samples (65536)
    const double delta_f = Fs / N; // frequency resolution

    std::cout << "Spectrum analysis configuration:" << std::endl;
    std::cout << "Fs = " << Fs << " Hz" << std::endl;
    std::cout << "N = " << N << " samples" << std::endl;
    std::cout << "T_total = " << (N / Fs) << " s" << std::endl;
    std::cout << "Delta_f = " << std::fixed << std::setprecision(5) << delta_f << " Hz" << std::endl;

    // --- 2. Generate base time-domain signal ---
    std::cout << "\nGenerating base signal..." << std::endl;
    std::vector<Complex> signal_base(N);
    double t;
    for (int i = 0; i < N; ++i)
    {
        t = (double)i / Fs; // current time (s)
        double s_t = 0.8 * std::sin(2 * PI * 103 * t) +
                     1.0 * std::sin(2 * PI * 107 * t) +
                     0.1 * std::sin(2 * PI * 115 * t);
        signal_base[i] = Complex(s_t, 0.0);
    }

    // --- 3. Process RECTANGULAR Window ---
    std::cout << "Processing Rectangular Window FFT..." << std::endl;
    std::vector<Complex> signal_rect = signal_base; // Create copy

    // No window function applied (equivalent to rectangular)
    fft(signal_rect); // In-place FFT

    // Open output file (truncate to clear old data)
    std::string rect_filename = get_filename("rect_spectrum", N);
    std::ofstream rect_outfile;
    rect_outfile.open(rect_filename, std::ios_base::trunc);
    if (!rect_outfile.is_open())
    {
        std::cerr << "Error: Could not open file " << rect_filename << std::endl;
        return 1;
    }

    std::cout << "Writing " << rect_filename << "..." << std::endl;
    for (int k = 0; k <= N / 2; ++k)
    {
        double freq = k * delta_f;
        double magnitude_raw = std::abs(signal_rect[k]);

        // ** CORRECT Rectangular window scaling **
        // DC (k=0): gain is N
        // AC (k>0): gain is N/2 for a sinusoid
        double amplitude_scaled = (k == 0) ? (magnitude_raw / N) : (magnitude_raw * 2.0 / N);

        rect_outfile << freq << "," << amplitude_scaled << "\n";
    }
    rect_outfile.close();

    // --- 4. Process HANN Window ---
    std::cout << "Processing Hann Window FFT..." << std::endl;
    std::vector<Complex> signal_hann = signal_base; // Create fresh copy

    apply_hann_window(signal_hann); // Apply Hann window
    fft(signal_hann);               // In-place FFT

    // Open output file (truncate to clear old data)
    std::string hann_filename = get_filename("hann_spectrum", N);
    std::ofstream hann_outfile;
    hann_outfile.open(hann_filename, std::ios_base::trunc);
    if (!hann_outfile.is_open())
    {
        std::cerr << "Error: Could not open file " << hann_filename << std::endl;
        return 1;
    }

    std::cout << "Writing " << hann_filename << "..." << std::endl;
    std::cout << "\n--- Hann Spectrum (amplitude > 0.01) ---" << std::endl;
    std::cout << "-------------------------------------------------------------" << std::endl;
    std::cout << "Bin (k) | Freq (Hz) | Amplitude (scaled)" << std::endl;
    std::cout << "-------------------------------------------------------------" << std::endl;

    for (int k = 0; k <= N / 2; ++k)
    {
        double freq = k * delta_f;
        double magnitude_raw = std::abs(signal_hann[k]);

        // ** CORRECT Hann window scaling **
        // DC (k=0): coherent gain is N * 0.5
        // AC (k>0): coherent gain is N * 0.5 / 2
        double amplitude_scaled = (k == 0) ? (magnitude_raw * 2.0 / N) : (magnitude_raw * 4.0 / N);

        hann_outfile << freq << "," << amplitude_scaled << "\n";

        // Print only significant peaks for Hann (Rect is too noisy)
        if (amplitude_scaled > 0.01)
        {
            std::cout << std::setw(7) << k << " | "
                      << std::setw(9) << std::fixed << std::setprecision(3) << freq << " | "
                      << std::setw(16) << std::fixed << std::setprecision(4) << amplitude_scaled
                      << std::endl;
        }
    }
    hann_outfile.close();

    std::cout << "\nAll processing complete." << std::endl;

    return 0;
}
