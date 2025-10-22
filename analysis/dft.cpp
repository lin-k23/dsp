#include "transforms.h"
#include <cmath>

// Define PI
const double PI = std::acos(-1.0);

std::vector<Complex> dft(const std::vector<Complex> &x)
{
    int N = x.size();
    std::vector<Complex> X(N); // Output vector

    // Iterate over each output frequency bin 'k'
    for (int k = 0; k < N; ++k)
    {
        X[k] = Complex(0.0, 0.0); // Initialize sum for this bin

        // Sum over each input time sample 'n'
        for (int n = 0; n < N; ++n)
        {
            // Calculate the twiddle factor: W_N^(kn) = exp(-j * 2*PI * k * n / N)
            double angle = -2.0 * PI * k * n / N;
            Complex W = std::polar(1.0, angle);

            X[k] += x[n] * W;
        }
    }
    return X;
}