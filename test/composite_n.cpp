#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

// Define aliases for clarity
using Complex = std::complex<double>;
using CVector = std::vector<Complex>;

// Define PI
const double PI = std::acos(-1.0);

/**
 * @brief Performs a 2-point DFT (Radix-2 butterfly) in-place.
 * out[0] = in[0] + in[1]
 * out[1] = in[0] - in[1]
 */
void dft_2_point(Complex& a, Complex& b) {
    Complex temp = a;
    a = temp + b;
    b = temp - b;
}

/**
 * @brief Performs a 3-point DFT (naive implementation).
 * out[0] = in[0] + in[1] * W3^0 + in[2] * W3^0
 * out[1] = in[0] + in[1] * W3^1 + in[2] * W3^2
 * out[2] = in[0] + in[1] * W3^2 + in[2] * W3^4
 */
void dft_3_point(Complex& a0, Complex& a1, Complex& a2) {
    // Twiddle factors for N=3
    Complex W3_1 = std::exp(Complex(0, -2.0 * PI / 3.0)); // W_3^1
    Complex W3_2 = std::exp(Complex(0, -4.0 * PI / 3.0)); // W_3^2

    Complex in0 = a0;
    Complex in1 = a1;
    Complex in2 = a2;

    a0 = in0 + in1 + in2;
    a1 = in0 + in1 * W3_1 + in2 * W3_2;
    a2 = in0 + in1 * W3_2 + in2 * W3_1 * W3_1; // W_3^4 = W_3^1
}

/**
 * @brief Computes a 6-point FFT using the Cooley-Tukey algorithm
 * decomposing N=6 into r1=2 and r2=3.
 *
 * This function treats the 1D input as a 2x3 matrix:
 *
 * x[n] -> x[n1*r2 + n0] -> matrix M[n1, n0]
 *
 * n=0 -> M[0, 0] = x[0]
 * n=1 -> M[0, 1] = x[1]
 * n=2 -> M[0, 2] = x[2]
 * n=3 -> M[1, 0] = x[3]
 * n=4 -> M[1, 1] = x[4]
 * n=5 -> M[1, 2] = x[5]
 *
 * @param x Input vector of size 6.
 * @return Output (DFT) vector of size 6.
 */
CVector fft_6_point(const CVector& x) {
    if (x.size() != 6) {
        throw std::invalid_argument("Input vector must be of size 6.");
    }

    // 1. First Pass: 3 column-wise DFTs of size 2
    //    We compute DFTs on (x[0], x[3]), (x[1], x[4]), (x[2], x[5])
    
    CVector G(6); // This holds the intermediate result (G[k1, n0])
    
    // Column 0 (n0=0)
    Complex g00 = x[0]; // M[0, 0]
    Complex g10 = x[3]; // M[1, 0]
    dft_2_point(g00, g10);
    G[0] = g00; // G[0, 0]
    G[3] = g10; // G[1, 0] (storing in 1D array)

    // Column 1 (n0=1)
    Complex g01 = x[1]; // M[0, 1]
    Complex g11 = x[4]; // M[1, 1]
    dft_2_point(g01, g11);
    G[1] = g01; // G[0, 1]
    G[4] = g11; // G[1, 1]

    // Column 2 (n0=2)
    Complex g02 = x[2]; // M[0, 2]
    Complex g12 = x[5]; // M[1, 2]
    dft_2_point(g02, g12);
    G[2] = g02; // G[0, 2]
    G[5] = g12; // G[1, 2]

    // 2. Twiddle Factor Multiplication
    //    H[k1, n0] = G[k1, n0] * W_N^(n0 * k1)
    
    CVector H = G; // Copy G to H
    
    // Twiddle factors for N=6
    // We only need W_6^1 and W_6^2
    Complex W6_1 = std::exp(Complex(0, -2.0 * PI / 6.0)); // W_6^1
    Complex W6_2 = std::exp(Complex(0, -4.0 * PI / 6.0)); // W_6^2
    
    // H[0, n0] are multiplied by W_6^0 = 1 (so no change)
    // H[1, 0] is multiplied by W_6^0 = 1 (no change)
    
    // H[1, 1] = G[1, 1] * W_6^(1*1)
    H[4] = G[4] * W6_1; // H[k1=1, n0=1] is at index 1*3 + 1 = 4
    
    // H[1, 2] = G[1, 2] * W_6^(2*1)
    H[5] = G[5] * W6_2; // H[k1=1, n0=2] is at index 1*3 + 2 = 5
    
    
    // 3. Second Pass: 2 row-wise DFTs of size 3
    //    We compute DFTs on (H[0,0], H[0,1], H[0,2])
    //    and (H[1,0], H[1,1], H[1,2])
    
    CVector X_matrix(6); // Final result in 2x3 matrix form
    
    // Row 0 (k1=0)
    Complex x00 = H[0]; // H[0, 0]
    Complex x01 = H[1]; // H[0, 1]
    Complex x02 = H[2]; // H[0, 2]
    dft_3_point(x00, x01, x02);
    X_matrix[0] = x00; // X[0, 0]
    X_matrix[1] = x01; // X[0, 1]
    X_matrix[2] = x02; // X[0, 2]

    // Row 1 (k1=1)
    Complex x10 = H[3]; // H[1, 0]
    Complex x11 = H[4]; // H[1, 1]
    Complex x12 = H[5]; // H[1, 2]
    dft_3_point(x10, x11, x12);
    X_matrix[3] = x10; // X[1, 0]
    X_matrix[4] = x11; // X[1, 1]
    X_matrix[5] = x12; // X[1, 2]

    // 4. Final Reordering
    //    The output X[k1, k0] is "row-major".
    //    We need to map it to the 1D output X[k]
    //    using k = k1 + k0 * r1 (k = k1 + k0*2)
    
    // k=0 -> k0=0, k1=0 -> X_matrix[0] (index 0)
    // k=1 -> k0=0, k1=1 -> X_matrix[3] (index 3)
    // k=2 -> k0=1, k1=0 -> X_matrix[1] (index 1)
    // k=3 -> k0=1, k1=1 -> X_matrix[4] (index 4)
    // k=4 -> k0=2, k1=0 -> X_matrix[2] (index 2)
    // k=5 -> k0=2, k1=1 -> X_matrix[5] (index 5)

    CVector X_out(6);
    X_out[0] = X_matrix[0]; // k=0
    X_out[1] = X_matrix[3]; // k=1
    X_out[2] = X_matrix[1]; // k=2
    X_out[3] = X_matrix[4]; // k=3
    X_out[4] = X_matrix[2]; // k=4
    X_out[5] = X_matrix[5]; // k=5

    return X_out;
}

// --- Main function to test the N=6 FFT ---
int main() {
    // Simple test input: impulse at n=1
    // x = [0, 1, 0, 0, 0, 0]
    // DFT should be [1, W6^1, W6^2, W6^3, W6^4, W6^5]
    CVector x = {0, 1, 0, 0, 0, 0};

    std::cout << "--- Input Vector (x) ---" << std::endl;
    for(const auto& val : x) {
        std::cout << val << std::endl;
    }

    // Compute the 6-point FFT
    CVector X = fft_6_point(x);

    std::cout << "\n--- Output Vector (X) ---" << std::endl;
    for (int k = 0; k < 6; ++k) {
        std::cout << "X[" << k << "] = " << X[k] << std::endl;
    }
    
    std::cout << "\n--- Verification (W_6^k) ---" << std::endl;
    for (int k = 0; k < 6; ++k) {
        Complex W6_k = std::exp(Complex(0, -2.0 * PI * k / 6.0));
        std::cout << "W_6^" << k << " = " << W6_k << std::endl;
    }

    return 0;
}