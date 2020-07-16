#ifndef MATH_HPP_
#define MATH_HPP_

#include <type_traits>
#include <cassert>
#include <limits>
#include <optional>
#include <algorithm>


namespace math {

/*
 * Calculating large product of two numbers modulo n
 * using modular arithmetic principles.
 * Function returns: (a x b) mod n
 */
template<typename T = int>
T mod_mul(T a, T b, T n);


/*
 * Calculating large powers of base modulo n
 * using modular arithmetic principles.
 * Function returns: pow(base, large_number) mod n.
 */
template<typename T = int>
T mod_pow(T e, T n, T base = 2);


/*
 * Checking if the input number is prime using Miller-Rabin test.
 */
template<typename T = int>
bool is_prime(T a, int n = 20);


/*
 * Generate a prime number within the specified range.
 * There may not be a prime number in the given interval.
 */
template<typename T = int>
std::optional<T> generate_prime(T min = 0,
                                T max = std::numeric_limits<T>::max() / 2);


/*
 * Finding greatest common divisor using the Euclidean algorithm.
 * To reduce the number of operations, the algorithm uses modulo function.
 */
template<typename T = int>
T gcd(T a, T b);


/*
 * Checking if the numbers are relatively prime
 * using the greatest common divisor (gcd).
 */
template<typename T>
bool coprime(T a, T b);


/*
 * Modular multiplication inverse using extended Euclidean algorithm.
 * The inverse of modulo for given number may not exist.
 */
template<typename T>
std::optional<T> mod_inv(T a, T mod);


template<typename T>
T mod_mul(T a, T b, const T n)
{
    static_assert(std::is_integral<T>::value, "T must be integral");
    assert(a >= 0 && b >= 0 && n > 0);

    T result {0};
    // bit mask used to test the set bits
    T mask {1};

    while (mask) {
        if (b & mask) {
            result = (result + a) % n;
        }

        a = (a << 1) % n;
        mask <<= 1;
    }

    return result;
}


template<typename T>
T mod_pow(T e, T n, const T base)
{
    assert(e >= 0 && n >= 0 && base >= 0);

    T result {1};
    // the rest of the next modulo powers
    T p {base};
    // bit mask used to test the set bits
    T mask {1};

    while (mask) {
        if (e & mask) {
            result = mod_mul(result, p, n);
        }

        p = mod_mul(p, p, n);
        mask <<= 1;
    }

    return result;
}


template<typename T = int>
bool is_prime(T a, int n)
{
    static_assert(std::is_integral<T>::value, "T must be integral");
    assert(a >= 0 && n > 0);

    if (a == 0 || a == 1) {
        return false;
    }

    if (a == 2 || a == 3) {
        return true;
    }

    // exponent of power 2 in the divior p - 1
    T s {0};
    // power multiplier 2 in the divior p - 1
    T d {a - 1};

    // remove the divisors of 2 in p - 1 and count them in s
    while (!(d % n)) {
        ++s;
        d /= 2;
    }

    // perform n Miller-Rabin tests
    for (T base {}, x{}; n; --n) {
        base = utility::rand(2, a - 2);
        // expression of Miller-Rabin sequence
        x = mod_pow(d, a, base);

        if (x == 1 || x == a - 1) {
            continue;
        }

        for (T j {1}; j < s && x != a - 1; ++j) {
            x = mod_mul(x, x, a);

            if (x == 1) {
                return false;
            }
        }

        if (x != a - 1) {
            return false;
        }
    }

    return true;
}


template<typename T = int>
std::optional<T> generate_prime(T min, T max)
{
    static_assert(std::is_integral<T>::value, "T must be integral");

    // rand start number form which the prime number will be searched
    T x = utility::rand(min, max);

    // check only odd numbers within specified range
    for (x = x % 2 == 0 ? x + 1 : x; x <= max; x += 2) {
        if (is_prime(x)) {
            return x;
        }
    }

    return std::nullopt;
}


template<typename T = int>
T gcd(T a, T b)
{
    static_assert(std::is_integral<T>::value, "T must be integral");
    assert(a >= 0 && b >= 0);
    assert(a || b);

    if (!a || !b) {
        return a == 0 ? b : a;
    }

    while (b) {
        T temp {b};

        b = std::move(a) % std::move(b);
        a = std::move(temp);
    }

    return a;
}


template<typename T>
bool coprime(T a, T b)
{
    if (gcd(a, b) == 1) {
        return true;
    }

    return false;
}


template<typename T>
std::optional<T> mod_inv(T a, T mod)
{
    // coefficents of equations
    T u {1};
    T x {0};
    T w {a};
    T z {mod};
    // total quotient
    T q {};

    while (w) {
        if (w < z) {
            std::swap(u, x);
            std::swap(w, z);
        }

        q = w / z;
        u -= q * x;
        w -= q * z;
    }

    if (z == 1) {
        if (x < 0) {
            x += mod;
        }

        return x;
    }

    return std::nullopt;
}

};

#endif // MATH_HPP_
