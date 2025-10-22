#include <iostream>
#include <vector>
#include <algorithm> // For std::swap

// Use long long for all calculations to prevent intermediate overflow
using ll = long long;
using Poly = std::vector<ll>;

// A common "FFT-friendly" prime: 998244353 = 119 * 2^23 + 1
// This allows for efficient NTT.
const ll MOD = 998244353;
// G is a primitive root modulo MOD.
const ll G = 3; 

/**
 * @brief Modular exponentiation (base^exp % MOD).
 */
ll power(ll base, ll exp) {
    ll res = 1;
    base %= MOD;
    while (exp > 0) {
        if (exp % 2 == 1) res = (res * base) % MOD;
        base = (base * base) % MOD;
        exp >>= 1;
    }
    return res;
}

/**
 * @brief Modular inverse (n^-1 % MOD) using Fermat's Little Theorem.
 */
ll inv(ll n) {
    return power(n, MOD - 2);
}

/**
 * @brief Performs the Number Theoretic Transform (NTT).
 * @param a The polynomial (vector of coefficients) to transform.
 * @param invert True for Inverse NTT (INTT), false for forward NTT.
 */
void ntt(Poly& a, bool invert) {
    int n = a.size();

    // 1. Bit-Reversal Permutation
    // This reorders the array for in-place butterfly operations.
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            std::swap(a[i], a[j]);
        }
    }

    // 2. Butterfly Operations (Cooley-Tukey DIT)
    for (int len = 2; len <= n; len <<= 1) {
        ll wlen = power(G, (MOD - 1) / len);
        if (invert) {
            wlen = inv(wlen);
        }
        for (int i = 0; i < n; i += len) {
            ll w = 1;
            for (int j = 0; j < len / 2; j++) {
                ll u = a[i + j];
                ll v = (w * a[i + j + len / 2]) % MOD;
                a[i + j] = (u + v) % MOD;
                a[i + j + len / 2] = (u - v + MOD) % MOD; // Add MOD to handle negatives
                w = (w * wlen) % MOD;
            }
        }
    }

    // 3. Scale by 1/n if it's an inverse transform
    if (invert) {
        ll n_inv = inv(n);
        for (ll& x : a) {
            x = (x * n_inv) % MOD;
        }
    }
}

/**
 * @brief Multiplies two polynomials (convolves) using NTT.
 */
Poly multiply(Poly a, Poly b) {
    // 1. Find the smallest power of 2 size
    int n = 1;
    int target_size = a.size() + b.size() - 1;
    while (n < target_size) {
        n <<= 1;
    }

    // 2. Pad with zeros
    a.resize(n, 0);
    b.resize(n, 0);

    // 3. Transform
    ntt(a, false);
    ntt(b, false);

    // 4. Point-wise multiplication
    Poly c(n);
    for (int i = 0; i < n; i++) {
        c[i] = (a[i] * b[i]) % MOD;
    }

    // 5. Inverse transform
    ntt(c, true);
    
    // 6. Resize to the correct final degree
    c.resize(target_size);
    return c;
}

/**
 * @brief Computes polynomial A(x)^n using fast power (binary exponentiation).
 */
Poly polyPower(Poly a, int n) {
    // Start with the identity polynomial P(x) = 1
    Poly res = {1}; 
    
    while (n > 0) {
        if (n % 2 == 1) {
            res = multiply(res, a);
        }
        a = multiply(a, a);
        n >>= 1;
    }
    return res;
}


// --- Main function to solve the C(n, k) problem ---
int main() {
    int n;
    std::cout << "Enter the value of n: ";
    std::cin >> n;

    if (n < 0) {
        std::cout << "n must be non-negative." << std::endl;
        return 1;
    }

    // Our base polynomial is (1 + x)
    // The coefficients are {1, 1} for x^0 and x^1
    Poly a = {1, 1};

    // Compute (1 + x)^n
    Poly result = polyPower(a, n);

    // The coefficients of the result are the C(n, k) values
    std::cout << "--- Combinations C(n, k) mod " << MOD << " ---" << std::endl;
    for (int k = 0; k < result.size(); ++k) {
        // result[k] is the coefficient of x^k, which is C(n, k)
        std::cout << "C(" << n << ", " << k << ") = " << result[k] << std::endl;
    }

    return 0;
}