#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <limits>
#include <type_traits>
#include <utility>
#include <cassert>

#include "string_utility.hpp"


namespace math {

template<typename T,
         typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
class Matrix;


template<typename T>
class Matrix<T>
{
public:
    using Iterator = typename std::vector<std::vector<T>>::iterator;
    using ConstIterator = typename std::vector<std::vector<T>>::const_iterator;

    Matrix(size_t rows, size_t elements,
           const std::vector<std::string> &matrix = {});

    Matrix(const Matrix<T> &matrix);
    Matrix<T>& operator=(const Matrix<T> &matrix);
    Matrix(Matrix<T> &&matrix);
    Matrix<T>& operator=(Matrix<T> &&matrix);

    T& operator()(size_t row);
    const T& operator()(size_t row) const;

    T& operator()(size_t row, size_t element);
    const T& operator()(size_t row, size_t element) const;

    std::vector<T>& operator[](size_t row);
    const std::vector<T>& operator[](size_t row) const;

    template<typename V>
    friend std::ostream& operator<<(std::ostream &out, const Matrix<V> &matrix);

    template<typename V>
    friend bool operator==(const Matrix<V> &m1, const Matrix<V> &m2);

    template<typename V>
    friend Matrix<V> operator+(const Matrix<V> &m1, const Matrix<V> &m2);

    template<typename V>
    friend Matrix<V>& operator+=(Matrix<V> &m1, const Matrix<V> &m2);

    template<typename V>
    friend Matrix<V> operator-(const Matrix<V> &m1, const Matrix<V> &m2);

    template<typename V>
    friend Matrix<V>& operator-=(Matrix<V> &m1, const Matrix<V> &m2);

    template<typename V>
    friend Matrix<V> operator-(Matrix<V> &matrix);

    template<typename V>
    friend Matrix<V> operator*(const Matrix<V> &matrix, V scalar);

    template<typename V>
    friend Matrix<V> operator*(V scalar, const Matrix<V> &matrix);

    template<typename V>
    friend Matrix<V>& operator*=(Matrix<V> &m1, V scalar);

    template<typename V>
    friend Matrix<V> operator*(const Matrix<V> &m1, const Matrix<V> &m2);

    template<typename V>
    friend Matrix<V> operator*=(Matrix<V> &m1, const Matrix<V> &m2);

    void transpose();

    Iterator begin() noexcept { return matrix_.begin(); }
    Iterator end() noexcept { return matrix_.end(); }
    ConstIterator begin() const noexcept { return matrix_.cbegin(); }
    ConstIterator end() const noexcept { return matrix_.cend(); }

    size_t rows() const { return rows_; };
    size_t elements() const { return elements_; };

    T min_element(size_t row) const;
    T min_element() const;
    T max_element(size_t row) const;
    T max_element() const;

    void allocate(size_t rows, size_t elements);
    void initialize(size_t rows, size_t elements,
                    const std::vector<std::string> &matrix);
    void initialize(const Matrix<T> &matrix);

    void resize(size_t rows, size_t elements);

    void fill_row(size_t row_num, const std::string &row);
private:
    std::vector<std::vector<T>> matrix_ {};
    size_t rows_ {0};
    size_t elements_ {0};
};


template<typename T>
Matrix<T>::Matrix(const size_t rows, const size_t elements,
                  const std::vector<std::string> &matrix)
    : rows_(rows), elements_(elements)
{
    initialize(rows, elements, matrix);
}


template<typename T>
Matrix<T>::Matrix(const Matrix<T> &matrix)
{
    initialize(matrix);
}


template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T> &matrix)
{
    initialize(matrix);

    return *this;
}


template<typename T>
Matrix<T>::Matrix(Matrix<T> &&matrix)
{
    initialize(std::move(matrix));
}


template<typename T>
Matrix<T>& Matrix<T>::operator=(Matrix<T> &&matrix)
{
    initialize(std::move(matrix));

    return *this;
}


template<typename T>
T& Matrix<T>::operator()(size_t row)
{
    return matrix_.at(row);
}


template<typename T>
const T& Matrix<T>::operator()(size_t row) const
{
    return matrix_.at(row);
}


template<typename T>
const T& Matrix<T>::operator()(const size_t row, const size_t element) const
{
    return matrix_.at(row).at(element);
}


template<typename T>
T& Matrix<T>::operator()(const size_t row, const size_t element)
{
    return matrix_.at(row).at(element);
}


template<typename T>
std::vector<T>& Matrix<T>::operator[](size_t row)
{
    return matrix_.at(row);
}


template<typename T>
const std::vector<T>& Matrix<T>::operator[](size_t row) const
{
    return matrix_.at(row);
}


template<typename T>
std::ostream& operator<<(std::ostream &os, const Matrix<T> &matrix)
{
    auto min_element_str_size {utility::to_string(matrix.min_element()).size()};
    auto max_element_str_size {utility::to_string(matrix.max_element()).size()};

    auto width {min_element_str_size > max_element_str_size
                ? min_element_str_size : max_element_str_size};

    for (const auto &row : matrix) {
        for (auto val : row) {
            os << std::setw(width) << std::right << val << " ";
        }
        os << '\n';
    }

    return os;
}


template<typename T>
bool operator==(const Matrix<T> &m1, const Matrix<T> &m2)
{
    if (&m1 == &m2) {
        return true;
    }

    if ((m1.rows_ != m2.rows_) || (m1.elements_ != m2.elements_)) {
        return false;
    }

    for (size_t row {0}; row < m1.rows_; ++row) {
        if (m1.matrix_[row] != m2.matrix_[row]) {
            return false;
        }
    }

    return true;
}


template<typename T>
Matrix<T> operator+(const Matrix<T> &m1, const Matrix<T> &m2)
{
    assert((m1.rows_ == m2.rows_) && (m1.elements_ == m2.elements_));

    Matrix<T> m(m1.rows_, m1.elements_);

    for (size_t row {0}; row < m.rows_; ++row) {
        for (size_t element {0}; element < m.elements_; ++element) {
            m(row, element) = m1(row, element) + m2(row, element);
        }
    }

    return m;
}


template<typename T>
Matrix<T>& operator+=(Matrix<T> &m1, const Matrix<T> &m2)
{
    m1 = m1 + m2;

    return m1;
}


template<typename T>
Matrix<T> operator-(const Matrix<T> &m1, const Matrix<T> &m2)
{
    assert((m1.rows_ == m2.rows_) && (m1.elements_ == m2.elements_));

    Matrix<T> m(m1.rows_, m1.elements_);

    for (size_t row {0}; row < m.rows_; ++row) {
        for (size_t element {0}; element < m.elements_; ++element) {
            m(row, element) = m1(row, element) - m2(row, element);
        }
    }

    return m;
}


template<typename T>
Matrix<T>& operator-=(Matrix<T> &m1, const Matrix<T> &m2)
{
    m1 = m1 - m2;

    return m1;
}


template<typename T>
Matrix<T> operator-(Matrix<T> &matrix)
{
    for (auto &row : matrix.matrix_) {
        for (auto &element : row) {
            element = -element;
        }
    }

    return matrix;
}


template<typename T>
Matrix<T> operator*(const Matrix<T> &matrix, T scalar)
{
    Matrix<T> tmp_matrix(matrix.rows_, matrix.elements_);

    for (size_t row {0}; row < matrix.rows_; ++row) {
        for (size_t element {0}; element < matrix.elements_; ++element) {
            tmp_matrix(row, element) = matrix(row, element) * scalar;
        }
    }

    return tmp_matrix;
}


template<typename T>
Matrix<T> operator*(T scalar, const Matrix<T> &matrix)
{
    return operator*(matrix, scalar);
}


template<typename T>
Matrix<T>& operator*=(Matrix<T> &m1, T scalar)
{
    m1 = m1 * scalar;

    return m1;
}


template<typename V>
Matrix<V> operator*(const Matrix<V> &m1, const Matrix<V> &m2)
{
    assert(m1.elements_ == m2.rows_);

    Matrix<V> m(m1.rows_, m1.elements_);

    for (size_t i {0}; i < m1.rows_; ++i) {
        std::vector<V> row(m1.elements_, 0);

        for (size_t j {0}; j < m1.elements_; ++j) {
            for (size_t k {0}; k < m2.elements_; ++k) {
                row[k] += m1[i][j] * m2[j][k];
            }
        }

        m[i] = row;
    }

    return m;
}


template<typename V>
Matrix<V> operator*=(Matrix<V> &m1, const Matrix<V> &m2)
{
    m1 = m1 * m2;

    return m1;
}


template<typename T>
void Matrix<T>::transpose()
{
    auto tmp_matrix(std::move(matrix_));

    std::swap(rows_, elements_);
    allocate(rows_, elements_);

    for (size_t row {0}; row < rows_; ++row) {
        for (size_t element {0}; element < elements_; ++element) {
            matrix_[row][element] = tmp_matrix[element][row];
        }
    }
}


template<typename T>
T Matrix<T>::min_element(const size_t row) const
{
    return *std::min_element(matrix_.at(row).begin(), matrix_.at(row).end());
}


template<typename T>
T Matrix<T>::min_element() const
{
    auto cur_min_element {std::numeric_limits<T>::max()};

    for (size_t row {0}; row < matrix_.size(); ++row) {
        cur_min_element = min_element(row) < cur_min_element
                          ? min_element(row) : cur_min_element;
    }

    return cur_min_element;
}


template<typename T>
T Matrix<T>::max_element(const size_t row) const
{
    return *std::max_element(matrix_.at(row).begin(), matrix_.at(row).end());
}


template<typename T>
T Matrix<T>::max_element() const
{
    auto cur_max_element {std::numeric_limits<T>::min()};

    for (size_t row {0}; row < matrix_.size(); ++row) {
        cur_max_element = max_element(row) > cur_max_element
                          ? max_element(row) : cur_max_element;
    }

    return cur_max_element;
}


template<typename T>
void Matrix<T>::allocate(const size_t rows, const size_t elements)
{
    rows_ = rows;
    elements_ = elements;

    matrix_.resize(rows);

    for (auto &row : matrix_) {
        row.resize(elements);
    }
}


template<typename T>
void Matrix<T>::initialize(const size_t rows, const size_t elements,
                           const std::vector<std::string> &matrix)
{
    allocate(rows, elements);

    for (size_t row {0}; row < rows_; ++row) {
        std::vector<T> row_vec {};

        utility::from_string(matrix[row], ' ', row_vec);

        matrix_[row] = row_vec;
    }
}


template<typename T>
void Matrix<T>::initialize(const Matrix<T> &matrix)
{
    allocate(matrix.rows(), matrix.elements());

    for (size_t row {0}; row < rows_; ++row) {
        matrix_[row] = matrix[row];
    }
}


template<typename T>
void Matrix<T>::resize(size_t rows, size_t elements)
{
    rows_ = rows;
    elements_ = elements;

    matrix_.resize(rows);

    for (auto &row : matrix_) {
        row.resize(elements);
    }
}


template<typename T>
void Matrix<T>::fill_row(const size_t row_num, const std::string &row)
{
    assert(row_num < rows_);

    std::vector<T> row_vec {};

    utility::from_string(row, ' ', row_vec);

    matrix_[row_num] = row_vec;
}

} // namespace math

#endif // MATRIX_HPP_
