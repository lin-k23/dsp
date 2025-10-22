#include <iostream>
#include <vector>
#include <complex>
#include <cmath>     // For std::log2
#include <algorithm> // For std::swap

// Define type aliases for easier reading
using Complex = std::complex<double>;
using CVector = std::vector<Complex>;

/**
 * @brief Reverses the bits of an index.
 * This implements the (n2 n1 n0) -> (n0 n1 n2) logic.
 * @param index The original index (e.g., 1, which is 001).
 * @param bits The total number of bits (e.g., 3 for N=8).
 * @return The bit-reversed index (e.g., 4, which is 100).
 */
unsigned int reverseBits(unsigned int index, unsigned int bits) {
    unsigned int reversedIndex = 0;
    for (unsigned int i = 0; i < bits; ++i) {
        // Check if the i-th bit of the original index is 1
        if ((index >> i) & 1) {
            // If it is, set the corresponding (bits - 1 - i)-th bit
            // in the reversed index.
            reversedIndex |= (1 << (bits - 1 - i));
        }
    }
    return reversedIndex;
}

/**
 * @brief Reorders a vector in-place using bit reversal.
 * @param data The vector of complex samples. Its size MUST be a power of 2.
 */
void bitReverseReorder(CVector& data) {
    const unsigned int N = data.size();

    // N must be a power of 2
    // This check works for N > 0: (N & (N - 1)) == 0
    if (N == 0 || (N & (N - 1)) != 0) {
        std::cerr << "Error: bitReverseReorder size must be a power of 2." << std::endl;
        return;
    }

    // Get the number of bits, e.g., log2(8) = 3
    unsigned int bits = static_cast<unsigned int>(std::log2(N));

    // Iterate through all indices
    for (unsigned int i = 0; i < N; ++i) {
        // Get the bit-reversed index
        unsigned int j = reverseBits(i, bits);

        //
        // IMPORTANT: Only swap if j > i
        // This prevents swapping elements twice (e.g., swapping 1 with 4,
        // and then later swapping 4 with 1, which would undo the work).
        //
        if (j > i) {
            std::swap(data[i], data[j]);
        }
    }
}

/**
 * @brief Helper function to print the vector for demonstration.
 */
void printVector(const CVector& v) {
    std::cout << "Index (Bin) | Reversed (Bin) | Value" << std::endl;
    std::cout << "---------------------------------------" << std::endl;
    
    unsigned int N = v.size();
    unsigned int bits = (N > 0) ? static_cast<unsigned int>(std::log2(N)) : 0;

    for (unsigned int i = 0; i < N; ++i) {
        unsigned int j = reverseBits(i, bits);
        std::cout << "  " << i << "   (..." << ((i >> 2) & 1) << ((i >> 1) & 1) << (i & 1) << ")"
                  << " |    " << j << "   (..." << ((j >> 2) & 1) << ((j >> 1) & 1) << (j & 1) << ")"
                  << " | " << v[i].real() << std::endl;
    }
}


// --- Main function to demonstrate the reordering ---
int main() {
    // Let's use N=8 (3 bits, as in your n2,n1,n0 example)
    const int N = 8;
    CVector data(N);

    // Initialize data with its original index (0, 1, 2, ... 7)
    // so we can clearly see where each element moves.
    for (int i = 0; i < N; ++i) {
        data[i] = Complex(i, 0);
    }

    std::cout << "--- Original Order ---" << std::endl;
    // We pass a copy for this print to show the "before" state
    printVector(CVector(data.begin(), data.end()));

    // Perform the in-place reordering
    bitReverseReorder(data);

    std::cout << "\n--- Bit-Reversed Order ---" << std::endl;
    printVector(data);

    return 0;
}