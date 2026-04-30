#include "matrix.h"
#include <iostream>

Matrix::Matrix(size_t r, size_t c) : rows(r), cols(c), data(r * c) {
    if (r == 0 || c == 0) throw std::invalid_argument("Нулевой размер матрицы.");
}

Fraction& Matrix::operator()(size_t r, size_t c) { return data[r * cols + c]; }
const Fraction& Matrix::operator()(size_t r, size_t c) const { return data[r * cols + c]; }

void Matrix::print() const {
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            std::cout << (*this)(i, j) << "\t";
        }
        std::cout << "\n";
    }
}

Matrix Matrix::operator+(const Matrix& other) const {
    Matrix res(rows, cols);
    for (size_t i = 0; i < data.size(); ++i) res.data[i] = data[i] + other.data[i];
    return res;
}

Matrix Matrix::operator-(const Matrix& other) const {
    Matrix res(rows, cols);
    for (size_t i = 0; i < data.size(); ++i) res.data[i] = data[i] - other.data[i];
    return res;
}

std::vector<Fraction> Matrix::operator*(const std::vector<Fraction>& vec) const {
    std::vector<Fraction> res(rows);
    for (size_t i = 0; i < rows; ++i) {
        Fraction sum(0);
        for (size_t j = 0; j < cols; ++j) sum += (*this)(i, j) * vec[j];
        res[i] = std::move(sum);
    }
    return res;
}

Matrix Matrix::operator*(const Matrix& other) const {
    Matrix res(rows, other.cols);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t k = 0; k < cols; ++k) {
            const Fraction& r_ik = (*this)(i, k);
            if (r_ik.isZero()) continue;
            for (size_t j = 0; j < other.cols; ++j) res(i, j) += r_ik * other(k, j);
        }
    }
    return res;
}

Matrix Matrix::multiplyBlock(const Matrix& other, size_t B) const {
    Matrix res(rows, other.cols);
    for (size_t ii = 0; ii < rows; ii += B) {
        for (size_t kk = 0; kk < cols; kk += B) {
            for (size_t jj = 0; jj < other.cols; jj += B) {
                for (size_t i = ii; i < std::min(ii + B, rows); ++i) {
                    for (size_t k = kk; k < std::min(kk + B, cols); ++k) {
                        const Fraction& r_ik = (*this)(i, k);
                        if (r_ik.isZero()) continue;
                        for (size_t j = jj; j < std::min(jj + B, other.cols); ++j) {
                            res(i, j) += r_ik * other(k, j);
                        }
                    }
                }
            }
        }
    }
    return res;
}

size_t Matrix::nextPowerOfTwo(size_t n) const {
    size_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

Matrix Matrix::padToPowerOfTwo(size_t newSize) const {
    Matrix res(newSize, newSize);
    for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; j < cols; ++j)
            res(i, j) = (*this)(i, j);
    return res;
}

Matrix Matrix::unpad(size_t origRows, size_t origCols) const {
    Matrix res(origRows, origCols);
    for (size_t i = 0; i < origRows; ++i)
        for (size_t j = 0; j < origCols; ++j)
            res(i, j) = (*this)(i, j);
    return res;
}

void Matrix::split(Matrix& A11, Matrix& A12, Matrix& A21, Matrix& A22) const {
    size_t half = rows / 2;
    for (size_t i = 0; i < half; ++i) {
        for (size_t j = 0; j < half; ++j) {
            A11(i, j) = (*this)(i, j);
            A12(i, j) = (*this)(i, j + half);
            A21(i, j) = (*this)(i + half, j);
            A22(i, j) = (*this)(i + half, j + half);
        }
    }
}

void Matrix::join(const Matrix& C11, const Matrix& C12, const Matrix& C21, const Matrix& C22) {
    size_t half = rows / 2;
    for (size_t i = 0; i < half; ++i) {
        for (size_t j = 0; j < half; ++j) {
            (*this)(i, j) = C11(i, j);
            (*this)(i, j + half) = C12(i, j);
            (*this)(i + half, j) = C21(i, j);
            (*this)(i + half, j + half) = C22(i, j);
        }
    }
}

Matrix Matrix::multiplyStrassen(const Matrix& other) const {
    size_t maxSize = std::max({rows, cols, other.rows, other.cols});
    size_t n = nextPowerOfTwo(maxSize);
    if (n <= 64) return (*this) * other;

    Matrix A = padToPowerOfTwo(n);
    Matrix B = other.padToPowerOfTwo(n);
    size_t half = n / 2;

    Matrix A11(half, half), A12(half, half), A21(half, half), A22(half, half);
    Matrix B11(half, half), B12(half, half), B21(half, half), B22(half, half);
    A.split(A11, A12, A21, A22);
    B.split(B11, B12, B21, B22);

    Matrix P1 = A11.multiplyStrassen(B12 - B22);
    Matrix P2 = (A11 + A12).multiplyStrassen(B22);
    Matrix P3 = (A21 + A22).multiplyStrassen(B11);
    Matrix P4 = A22.multiplyStrassen(B21 - B11);
    Matrix P5 = (A11 + A22).multiplyStrassen(B11 + B22);
    Matrix P6 = (A12 - A22).multiplyStrassen(B21 + B22);
    Matrix P7 = (A11 - A21).multiplyStrassen(B11 + B12);

    Matrix C(n, n);
    C.join(P5 + P4 - P2 + P6, P1 + P2, P3 + P4, P5 + P1 - P3 - P7);
    return C.unpad(rows, other.cols);
}

void Matrix::swapRows(size_t r1, size_t r2) {
    if (r1 == r2) return;
    for (size_t j = 0; j < cols; ++j) std::swap((*this)(r1, j), (*this)(r2, j));
}

void Matrix::pivot(size_t col, size_t startRow, bool verbose) {
    size_t maxRow = startRow;
    Fraction zero(0);
    Fraction minusOne(-1);
    for (size_t i = startRow + 1; i < rows; ++i) {
        Fraction current = (*this)(i, col);
        Fraction maxElem = (*this)(maxRow, col);
        Fraction currentAbs = (current < zero) ? (current * minusOne) : current;
        Fraction maxAbs = (maxElem < zero) ? (maxElem * minusOne) : maxElem;
        if (currentAbs > maxAbs) maxRow = i;
    }
    if (startRow != maxRow) {
        if (verbose) std::cout << "  [Пивотинг] Строки " << startRow + 1 << " <-> " << maxRow + 1 << "\n";
        swapRows(startRow, maxRow);
    }
    if ((*this)(startRow, col).isZero()) throw std::runtime_error("Матрица вырождена.");
}

std::vector<Fraction> Matrix::solveGauss(const std::vector<Fraction>& b, bool verbose) const {
    Matrix aug(rows, cols + 1);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) aug(i, j) = (*this)(i, j);
        aug(i, cols) = b[i];
    }
    for (size_t i = 0; i < rows; ++i) {
        aug.pivot(i, i, verbose);
        for (size_t k = i + 1; k < rows; ++k) {
            Fraction factor = aug(k, i) / aug(i, i);
            for (size_t j = i; j <= cols; ++j) aug(k, j) -= factor * aug(i, j);
        }
    }
    std::vector<Fraction> x(rows);
    for (int i = rows - 1; i >= 0; --i) {
        x[i] = aug(i, cols);
        for (size_t j = i + 1; j < cols; ++j) x[i] -= aug(i, j) * x[j];
        x[i] /= aug(i, i);
    }
    return x;
}

std::vector<Fraction> Matrix::solveJordanGauss(const std::vector<Fraction>& b, bool verbose) const {
    Matrix aug(rows, cols + 1);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) aug(i, j) = (*this)(i, j);
        aug(i, cols) = b[i];
    }
    for (size_t i = 0; i < rows; ++i) {
        aug.pivot(i, i, verbose);
        Fraction diag = aug(i, i);
        for (size_t j = i; j <= cols; ++j) aug(i, j) /= diag;
        for (size_t k = 0; k < rows; ++k) {
            if (k != i) {
                Fraction factor = aug(k, i);
                for (size_t j = i; j <= cols; ++j) aug(k, j) -= factor * aug(i, j);
            }
        }
    }
    std::vector<Fraction> x(rows);
    for (size_t i = 0; i < rows; ++i) x[i] = aug(i, cols);
    return x;
}

Matrix Matrix::inverseGauss(bool verbose) const {
    if (rows != cols) throw std::invalid_argument("Матрица должна быть квадратной.");
    Matrix aug(rows, cols * 2);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) aug(i, j) = (*this)(i, j);
        aug(i, cols + i) = Fraction(1);
    }
    for (size_t i = 0; i < rows; ++i) {
        aug.pivot(i, i, verbose);
        for (size_t k = i + 1; k < rows; ++k) {
            Fraction factor = aug(k, i) / aug(i, i);
            for (size_t j = i; j < cols * 2; ++j) aug(k, j) -= factor * aug(i, j);
        }
    }
    Matrix inv(rows, cols);
    for (size_t c = 0; c < cols; ++c) {
        for (int i = rows - 1; i >= 0; --i) {
            Fraction sum = aug(i, cols + c);
            for (size_t j = i + 1; j < cols; ++j) sum -= aug(i, j) * inv(j, c);
            inv(i, c) = sum / aug(i, i);
        }
    }
    return inv;
}

Matrix Matrix::inverseJordanGauss(bool verbose) const {
    if (rows != cols) throw std::invalid_argument("Матрица должна быть квадратной.");
    Matrix aug(rows, cols * 2);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) aug(i, j) = (*this)(i, j);
        aug(i, cols + i) = Fraction(1);
    }
    for (size_t i = 0; i < rows; ++i) {
        aug.pivot(i, i, verbose);
        Fraction diag = aug(i, i);
        for (size_t j = i; j < cols * 2; ++j) aug(i, j) /= diag;
        for (size_t k = 0; k < rows; ++k) {
            if (k != i) {
                Fraction factor = aug(k, i);
                for (size_t j = i; j < cols * 2; ++j) aug(k, j) -= factor * aug(i, j);
            }
        }
    }
    Matrix inv(rows, cols);
    for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; j < cols; ++j) inv(i, j) = aug(i, cols + j);
        return inv;
}
