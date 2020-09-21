# Math Libray

A set of mathematical objects and functions for performing
mathematical operations in C++.

---
## Description

The library consists of 2 modules:
* math operations (math_opeations.hpp header)
* matrix object (matrix.hpp header)

All of them provide specific functions and methods described below.

### Math operations module

This module contains a collection of template methods for
performing mathmetical operations.

* <strong>T mod_mul(T a, T b, T n)</strong> - calculating large product of two numbers modulo n
                     using modular arithmetic principles.
* <strong>T mod_pow(T e, T n, T base = 2)</strong> - calculating large powers of base modulo n
                            using modular arithmetic principles.
* <strong>bool is_prime(T a, int n = 20)</strong> - checking if the input number is
                                   prime using Miller-Rabin test.
* <strong>std::optional&lt;T&gt; generate_prime(T min = 0, T max = std::numeric_limits<T>::max() / 2)</strong> - 
generate a prime number within the specified range.
There may not be a prime number in the given interval.
* <strong>T gcd(T a, T b)</strong> - finding greatest common divisor using the Euclidean algorithm.
                    To reduce the number of operations, the algorithm uses modulo function.
* <strong>bool coprime(T a, T b)</strong> - checking if the numbers are relatively prime
                           using the greatest common divisor (gcd).
* <strong>std::optional&lt;T&gt; mod_inv(T a, T mod)</strong> - modular multiplication inverse using
                                         extended Euclidean algorithm.
                                         The inverse of modulo for given number may not exist.

### Matrix module
This module contains a template Matrix class for performing matrix operations.

Constructors:
* <strong>Matrix()</strong> 
* <strong>Matrix(size_t size, T init_value = 0)</strong>
* <strong>Matrix(size_t rows, size_t elements, T init_value = 0)</strong>
* <strong>Matrix(size_t size, const std::vector<std::string> &matrix)</strong> 
* <strong>Matrix(size_t size, std::vector<std::string> &&matrix)</strong>
* <strong>Matrix(size_t rows, size_t elements, const std::vector<std::string> &matrix)</strong>
* <strong>Matrix(size_t rows, size_t elements, std::vector<std::string> &&matrix)</strong>
* <strong>suport for copy construtor & move constructor & assigment operator & move assigment operator</strong>

Operators:
* <strong>T& operator()(size_t row) const</strong>
* <strong>T& operator()(size_t row, size_t element) const</strong>
* <strong>std::vector<T>& operator[](size_t row) const</strong>
* <strong>friend std::ostream& operator<<(std::ostream &out, const Matrix<V> &matrix)</strong>
* <strong>friend bool operator==(const Matrix<V> &m1, const Matrix<V> &m2)</strong>
* <strong>friend Matrix<V> operator+(const Matrix<V> &m1, const Matrix<V> &m2)</strong> 
* <strong>friend Matrix<V>& operator+=(Matrix<V> &m1, const Matrix<V> &m)</strong> 
* <strong>friend Matrix<V>& operator-(Matrix<V> &m1, const Matrix<V> &m)</strong>
* <strong>friend Matrix<V>& operator-=(Matrix<V> &m1, const Matrix<V> &m)</strong>
* <strong>friend Matrix<V> operator*(const Matrix<V> &matrix, V scalar</strong> 
* <strong>friend Matrix<V> operator*(V scalar, const Matrix<V> &matrix)</strong>
* <strong>friend Matrix<V>& operator*=(Matrix<V> &m1, V scalar)</strong>
* <strong>friend Matrix<V> operator*(const Matrix<V> &m1, const Matrix<V> &m2)</strong>
* <strong>friend Matrix<V> operator*=(Matrix<V> &m1, const Matrix<V> &m2)</strong>
  
Methods:
* <strong>void transpose()</strong>
* <strong>support for foreach iteration</strong>
* <strong>size_t rows() const</strong>
* <strong>size_t elements() const</strong>
* <strong>std::optional<T> min_element(size_t row) const</strong>
* <strong>std::optional<T> max_element(size_t row) const</strong>
* <strong>void resize(size_t size)</strong>
* <strong>void resize(size_t rows, size_t elements)</strong>
* <strong>void initialize(size_t size, V &&matrix)</strong>
* <strong>void initialize(size_t rows, size_t elements, V &&matrix)</strong>
* <strong>void initialize(const Matrix<T> &matrix)</strong>
* <strong>void initialize(Matrix<T> &&matrix)</strong>
* <strong>void fill_row(size_t row_num, V &&row)</strong>
---
## Technology
* C++17
* cmake
* make
---

## Requirements
* Operation system: Windows, Linux, macOS
* C++17 compiler
* cmake tool installed
* make tool installed
---

## External library
* utility-library: https://github.com/dmarcini/utility-library.git
---

## Building & Installing
Example for Linux system.

* install external <em>utility</em> library (like below)
* install math library:
```
git clone library-url
cd path-to-clone-library
mkdir build
cd build
cmake ..
sudo make install
```

## Usage
To use the library you have to include appropriate headers e.g.:
* #include "matrix.hpp" to use Matrix class
* #include "math_operations.hpp" to use math operations

All classes and funtions are placed in <em>math</em> namespace.
<br>
For example to create a Matrix object you have to type:
```
math::Matrix matrix(args);
```
