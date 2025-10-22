#include "transforms.h"
#include <cmath>

// Define PI
const double PI = std::acos(-1.0);

/**
 * @brief Apply a Hann (Hanning) window to the signal in-place
 * w[n] = 0.5 * (1 - cos(2*PI*n / (N-1)))
 */
void apply_hann_window(std::vector<Complex> &x)
{
    int N = x.size();
    if (N == 0)
        return;

    for (int n = 0; n < N; ++n)
    {
        double window_multiplier = 0.5 * (1.0 - std::cos(2.0 * PI * n / (N - 1)));

        // Multiply whole complex sample by the window (works for real input stored in complex)
        x[n] *= window_multiplier;
    }
}
