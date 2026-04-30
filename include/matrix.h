#pragma once
#include <vector>
#include <stdexcept>
#include <algorithm>
#include "fraction.h"
class Matrix {
private:
    size_t rows;
    size_t cols;
    std::vector<Fraction> data;
    size_t nextPowerOfTwo(size_t n) const;
    Matrix padToPowerOfTwo(size_t newSize) const;
    Matrix unpad(size_t origRows, size_t origCols) const;
    void split(Matrix& A11, Matrix& A12, Matrix& A21, Matrix& A22) const;
    void join(const Matrix& C11, const Matrix& C12, const Matrix& C21, const Matrix& C22);
    void swapRows(size_t r1, size_t r2);
    void pivot(size_t col, size_t startRow, bool verbose = false);
public:
    Matrix(size_t r, size_t c);
    Fraction& operator()(size_t r, size_t c);
    const Fraction& operator()(size_t r, size_t c) const;
    size_t getRows() const { return rows; }
    size_t getCols() const { return cols; }
    void print() const;
    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    std::vector<Fraction> operator*(const std::vector<Fraction>& vec) const;
    Matrix operator*(const Matrix& other) const;
    Matrix multiplyBlock(const Matrix& other, size_t blockSize = 32) const;
    Matrix multiplyStrassen(const Matrix& other) const;
    // test1
    std::vector<Fraction> solveGauss(const std::vector<Fraction>& b, bool verbose = false) const;
    std::vector<Fraction> solveJordanGauss(const std::vector<Fraction>& b, bool verbose = false) const;
    Matrix inverseGauss(bool verbose = false) const;
    Matrix inverseJordanGauss(bool verbose = false) const;
};
