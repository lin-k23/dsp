#include <iostream>
#include <vector>
#include <numeric>   // For std::accumulate
#include <algorithm> // For std::swap, std::reverse
#include <complex>
#include <iomanip>   // For std::setw
#include <stdexcept> // For std::runtime_error

/**
 * @brief Calculates the digit-reversed index for a mixed-radix DIT FFT.
 *
 * This function maps an input index 'n' to its permuted index 'j'.
 * It implements the DIT mapping:
 * Decompose n = n0 + n1*r2 + n2*r2*r3 + ...
 * Recompose j = n_m + ... + n2*r1*r2 + n1*r1 + n0
 *
 * A more direct way:
 * Decompose n using {rm, r(m-1), ..., r1}
 * Recompose j using {r1, r2, ..., rm}
 *
 * @param n The input index (0 to N-1).
 * @param radices The vector of factors {r1, r2, ..., rm}.
 * @return The permuted (digit-reversed) index j.
 */
int getDigitReversedIndex(int n, const std::vector<int> &radices)
{
    int m = radices.size();
    std::vector<int> digits(m);
    int temp_n = n;

    // 1. Decompose n using reversed radices {rm, r(m-1), ..., r1}
    for (int i = 0; i < m; ++i)
    {
        int r = radices[m - 1 - i]; // r_m, r_{m-1}, ...
        digits[i] = temp_n % r;
        temp_n /= r;
    }

    // 2. Recompose j using forward radices {r1, r2, ..., rm}
    int j = 0;
    for (int i = 0; i < m; ++i)
    {
        j = j * radices[i] + digits[i];
    }
    return j;
}

/**
 * @brief Generates the full permutation table for a given N and radices.
 * @param N The total size of the FFT.
 * @param radices The vector of factors.
 * @return A std::vector<int> where table[n] = j.
 */
std::vector<int> generatePermutationTable(int N, const std::vector<int> &radices)
{
    std::vector<int> table(N);
    for (int n = 0; n < N; ++n)
    {
        table[n] = getDigitReversedIndex(n, radices);
    }
    return table;
}

/**
 * @brief Reorders a vector of data in-place using the permutation table.
 *
 * This function templated to work with any data type (e.g., complex, int).
 * It uses the efficient swap-if-greater logic, which works because the
 * DIT digit-reversal map is an involution (f(f(n)) = n).
 *
 * @param data The std::vector of data to be reordered.
 * @param table The permutation table from generatePermutationTable.
 */
template <typename T>
void reorderInPlace(std::vector<T> &data, const std::vector<int> &table)
{
    if (data.size() != table.size())
    {
        throw std::runtime_error("Data and table sizes do not match.");
    }

    int N = data.size();
    for (int n = 0; n < N; ++n)
    {
        int j = table[n];
        // Only swap if j > n to ensure each pair is swapped exactly once
        if (j > n)
        {
            std::swap(data[n], data[j]);
        }
    }
}

// --- Testbench ---

/**
 * @brief Runs a single test case for the digit-reversal algorithm.
 */
void runTest(const std::string &testName, const std::vector<int> &radices)
{
    std::cout << "---------------------------------------\n";
    std::cout << "Test Case: " << testName << "\n";
    std::cout << "Radices: { ";
    for (int r : radices)
        std::cout << r << " ";
    std::cout << "}\n";

    // Calculate N
    int N = 1;
    for (int r : radices)
    {
        N *= r;
    }

    // 1. Generate permutation table
    std::vector<int> table = generatePermutationTable(N, radices);

    std::cout << "\nPermutation Map (n -> j):\n";
    for (int n = 0; n < N; ++n)
    {
        std::cout << "  " << std::setw(2) << n << " -> " << std::setw(2) << table[n] << "\n";
    }

    // 2. Create test data (0, 1, 2, ..., N-1)
    std::vector<int> data(N);
    std::iota(data.begin(), data.end(), 0); // Fills data with 0, 1, 2...

    std::cout << "\nOriginal Data:\n  ";
    for (int val : data)
        std::cout << val << " ";
    std::cout << "\n";

    // 3. Perform in-place reordering
    reorderInPlace(data, table);

    std::cout << "\nReordered Data (in-place):\n  ";
    for (int val : data)
        std::cout << val << " ";
    std::cout << "\n";

    // 4. Verification
    bool passed = true;
    for (int n = 0; n < N; ++n)
    {
        int j = table[n];
        // The data that was originally at index 'n' should now be at index 'j'.
        // After swapping, the data at index 'j' should be 'n'.
        if (data[j] != n)
        {
            passed = false;
        }
    }

    std::cout << "\nVerification (data[j] == n): " << (passed ? "PASSED" : "FAILED") << "\n";
    std::cout << "---------------------------------------\n";
}

int main()
{
    // // Test 1: Radix-2 (N=8) -> Standard "bit-reversal"
    // runTest("Radix-2 (N=8)", {2, 2, 2});

    // // Test 2: Mixed-Radix (N=6) -> The {2, 3} case
    // runTest("Mixed-Radix (N=6)", {2, 3});

    // // Test 3: Mixed-Radix (N=12) -> {2, 2, 3} case
    // runTest("Mixed-Radix (N=12)", {2, 2, 3});

    runTest("Mixed-Radix (N=9)", {3, 3});

    return 0;
}