/// @file Matrix.cpp
/// @brief Source file for Matrix class
/// @date 2025-04-09
/// @author Martín Hernández Tonzán

#include "../include/Matrix.h"
#include <iostream>
#include <iomanip>
#include <cmath>

Matrix::Matrix(int fil, int col) : fil(fil), col(col) {
    initMatrix();
}

Matrix::Matrix(int fil, int col, double v[], int n) : fil(fil), col(col) {
    initMatrix();

    int k = 0;

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++) {
            if (k < n)
                matrix[i][j] = v[k++];
            else
                matrix[i][j] = 0;
        }
}

Matrix::Matrix(const Matrix &m) : fil(m.fil), col(m.col) {
    initMatrix();

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            matrix[i][j] = m.matrix[i][j];
}

Matrix::~Matrix() {
    for (int i = 0; i < fil; i++)
        delete[] matrix[i];

    delete[] matrix;
}

void Matrix::initMatrix() {
    matrix = new double *[fil];
    for (int i = 0; i < fil; i++)
        matrix[i] = new double[col];

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            matrix[i][j] = 0.0;
}

int Matrix::nfils() const {
    return fil;
}

int Matrix::ncols() const {
    return col;
}

Matrix &Matrix::operator=(const Matrix &matrix2) {
    if (this != &matrix2) {
        // If dimensions do not match, reserve memory for the missing rows/columns
        if (fil != matrix2.fil || col != matrix2.col) {
            for (int i = 0; i < fil; ++i) delete[] matrix[i];
            delete[] matrix;
            fil = matrix2.fil;
            col = matrix2.col;
            initMatrix();
        }

        for (int i = 0; i < fil; i++)
            for (int j = 0; j < col; j++)
                matrix[i][j] = matrix2.matrix[i][j];
    }

    return *this;
}

Matrix Matrix::operator+(const Matrix &matrix2) {
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] + matrix2.matrix[i][j];

    return result;
}

Matrix Matrix::operator-(const Matrix &matrix2) {
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] - matrix2.matrix[i][j];

    return result;
}

Matrix Matrix::operator*(const Matrix &matrix2) {
    if (this->col != matrix2.fil) {
        printf("\nfailure in %s() line %d\n", __func__, __LINE__);
        throw std::invalid_argument("Incompatible matrix dimensions for multiplication.");
    }
    Matrix result(fil, matrix2.col);

    for (int i = 0; i < this->fil; i++) {
        for (int j = 0; j < matrix2.col; j++) {
            for (int k = 0; k < this->col; k++) {
                result.matrix[i][j] = result.matrix[i][j] + this->matrix[i][k] * matrix2.matrix[k][j];
            }
        }
    }

    return result;
}

Matrix Matrix::operator-() const {
    Matrix result(fil, col);
    for (int i = 1; i <= fil; ++i)
        for (int j = 1; j <= col; ++j)
            result(i, j) = -this->operator()(i, j);
    return result;
}

double &Matrix::operator()(const int i, const int j) const {
    return matrix[i - 1][j - 1];
}

double &Matrix::operator()(const int i) const {
    if (fil == 1 && i >= 1 && i <= col) {
        // Row vector
        return matrix[0][i - 1];
    } else if (col == 1 && i >= 1 && i <= fil) {
        // Column vector
        return matrix[i - 1][0];
    } else {
        std::cerr << "Access with (int) on matrix of size " << fil << " x " << col
                  << " and index i = " << i << std::endl;
        throw std::out_of_range("Invalid index (neither row or column).");
    }
}

// --> Matrix scalar operations

Matrix Matrix::operator*(double scalar) const {
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] * scalar;

    return result;
}

Matrix Matrix::operator/(double scalar) const {
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] / scalar;

    return result;
}

Matrix Matrix::operator+(double scalar) const {
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] + scalar;

    return result;
}
Matrix Matrix::operator-(double scalar) const {
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] - scalar;

    return result;
}

Matrix operator*(double scalar, const Matrix &m) {
    return m * scalar;
}

Matrix operator/(double scalar, const Matrix &m) {
    return m / scalar;
}

Matrix operator+(double scalar, const Matrix &m) {
    return m + scalar;
}
Matrix operator-(double scalar, const Matrix &m) {
    return m - scalar;
}

Matrix Matrix::transpose() const {
    Matrix temp(col, fil);
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            temp.matrix[j][i] = matrix[i][j];
    return temp;
}

Matrix Matrix::inv() {
    if (fil != col)
        throw std::runtime_error("Matrix must be square to invert.");

    int n = fil;
    Matrix result(n, n);
    Matrix temp(*this);

    for (int i = 0; i < n; i++)
        result.matrix[i][i] = 1.0;

    for (int i = 0; i < n; i++) {
        double pivot = temp.matrix[i][i];
        if (fabs(pivot) < 1e-12)
            throw std::runtime_error("Matrix is singular or nearly singular.");

        for (int j = 0; j < n; j++) {
            temp.matrix[i][j] /= pivot;
            result.matrix[i][j] /= pivot;
        }

        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = temp.matrix[k][i];
                for (int j = 0; j < n; j++) {
                    temp.matrix[k][j] -= factor * temp.matrix[i][j];
                    result.matrix[k][j] -= factor * result.matrix[i][j];
                }
            }
        }
    }

    return result;
}

Matrix Matrix::eye(int n) {
    Matrix I(n, n);
    for (int i = 0; i < n; ++i)
        I.matrix[i][i] = 1.0;
    return I;
}

Matrix Matrix::subvector(int from, int to) const {
    if (from < 1 || to > (fil == 1 ? col : fil) || from > to)
        throw std::out_of_range("Invalid sub vector range.");

    int size = to - from + 1;
    Matrix result(1, size);

    for (int i = 0; i < size; ++i)
        result(1, i + 1) = (*this)(from + i);

    return result;
}

Matrix Matrix::concatRow(const Matrix &other) const {
    if (this->fil != other.fil)
        throw std::invalid_argument("Matrices must have the same number of rows to concatenate horizontally.");

    Matrix result(this->fil, this->col + other.col);

    for (int i = 0; i < this->fil; ++i) {
        for (int j = 0; j < this->col; ++j)
            result.matrix[i][j] = this->matrix[i][j];
        for (int j = 0; j < other.col; ++j)
            result.matrix[i][this->col + j] = other.matrix[i][j];
    }

    return result;
}

Matrix Matrix::getColumn(int column) const {
    Matrix result(fil, 1);
    for (int i = 1; i <= fil; ++i)
        result.matrix[i - 1][0] = this->matrix[i - 1][column - 1];
    return result;
}

double Matrix::norm(const Matrix &m) {
    double sum = 0.0;
    for (int i = 0; i < m.fil; ++i) {
        for (int j = 0; j < m.col; ++j) {
            sum += m.matrix[i][j] * m.matrix[i][j];
        }
    }
    return sqrt(sum);
}

double Matrix::dot(const Matrix &m1, const Matrix &m2) {
    if (m1.fil != m2.fil || m1.col != m2.col) {
        throw std::invalid_argument("Matrix dimensions must match for dot product.");
    }

    double sum = 0.0;
    for (int i = 0; i < m1.fil; ++i) {
        for (int j = 0; j < m1.col; ++j) {
            sum += m1.matrix[i][j] * m2.matrix[i][j];
        }
    }
    return sum;
}

Matrix Matrix::cross(const Matrix &m1, const Matrix &m2) {
    if (m1.fil != 3 || m2.fil != 3 || m1.col != 1 || m2.col != 1) {
        throw std::invalid_argument("Matrix is not 3x1.");
    }

    Matrix result(3, 1);

    result(1) = m1(2) * m2(3) - m1(3) * m2(2);
    result(2) = m1(3) * m2(1) - m1(1) * m2(3);
    result(3) = m1(1) * m2(2) - m1(2) * m2(1);

    return result;
}

Matrix Matrix::range(int from, int to, int step) {
    int size = ((to - from) / step) + 1;
    Matrix result(1, size);
    int j = 1;
    for (int i = from; i <= to; i += step) {
        result(1, j) = i;
        j++;
    }
    return result;
}

void Matrix::print() {
    for (int i = 0; i < fil; i++) {
        for (int j = 0; j < col; j++) {
            std::cout << std::fixed << std::setprecision(14) << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}



