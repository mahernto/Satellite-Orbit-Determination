/// @file Matrix.h
/// @brief Header file for Matrix class
/// @date 2025-04-09
/// @author Martín Hernández Tonzán

#ifndef _MATRIX_
#define _MATRIX_

class Matrix
{
    public:
        /// @brief Constructs a Matrix initialized with zeros
        /// @param[in] fil number of rows
        /// @param[in] col number of columns
        Matrix(int fil, int col);

        /// @brief Constructs a Matrix using a double array
        /// @param[in] fil number of rows
        /// @param[in] col number of columns
        /// @param[in] v double array which stores the matrix elements
        /// @param[in] n size of the array v
        Matrix(int fil, int col, double v[], int n);

        /// @brief Copy constructor
        /// @param[in] m Matrix to be copied
        Matrix(const Matrix& m);

        /// @brief Destructor that releases allocated memory.
        ~Matrix();

        /// @brief rows of the Matrix getter
        /// @return number of rows of the Matrix
        int nfils() const;
        /// @brief columns of the Matrix getter
        /// @return number of columns of the Matrix
        int ncols() const;

        /// @brief Overloads Matrix assignment operator.
        /// @param[in] matrix2 Matrix to assign from.
        /// @return reference to the assigned matrix.
        Matrix& operator=(const Matrix& matrix2);

        /// @brief Adds two Matrix objects
        /// @param[in] matrix2 Matrix to add with
        /// @return sum of the both Matrix
        Matrix  operator+(const Matrix& matrix2);

        /// @brief Subtracts one Matrix for another
        /// @param[in] matrix2 Matrix to subtract with
        /// @return Matrix after subtraction
        Matrix  operator-(const Matrix& matrix2);

        /// @brief Multiplies two Matrix objects
        /// @param[in] matrix2 Matrix to multiply with
        /// @return product of the both Matrix
        Matrix  operator*(const Matrix& matrix2);

        /// @brief Opposite of the matrix
        /// @return opposite matrix
        Matrix  operator-() const;

        /// @brief Accesses an element of the matrix
        /// @param[in] i row index
        /// @param[in] j column index
        /// @return reference to the element at position (i, j)
        double& operator()(const int i, const int j) const;

        /// @brief Accesses an element of the first row or columns of the matrix
        /// @param[in] j column or row index
        /// @return reference to the element at position (j)
        double& operator()(const int j) const;

        // Matrix scalar operations

        /// @brief Multiplies a Matrix and a scalar
        /// @param[in] scalar double value to multiply with
        /// @return product of Matrix and scalar
        Matrix operator*(double scalar) const;
        /// @brief Divides a Matrix and a scalar
        /// @param[in] scalar double value to divide with
        /// @return result of the division of Matrix and scalar
        Matrix operator/(double scalar) const;
        /// @brief Adds a Matrix and a scalar
        /// @param[in] scalar double value to add with
        /// @return addition of Matrix and scalar
        Matrix operator+(double scalar) const;
        /// @brief Subtracts a Matrix and a scalar
        /// @param[in] scalar double value to subtract with
        /// @return subtraction of Matrix and scalar
        Matrix operator-(double scalar) const;
        /// @brief Multiplies a Matrix and a scalar (scalar on the left)
        /// @param[in] scalar double value to multiply with
        /// @return product of Matrix and scalar
        friend Matrix operator*(double scalar, const Matrix& m);
        /// @brief Divides a Matrix and a scalar (scalar on the left)
        /// @param[in] scalar double value to divide with
        /// @return result of the division of Matrix and scalar
        friend Matrix operator/(double scalar, const Matrix& m);
        /// @brief Adds a Matrix and a scalar (scalar on the left)
        /// @param[in] scalar double value to add with
        /// @return addition of Matrix and scalar
        friend Matrix operator+(double scalar, const Matrix& m);
        /// @brief Subtracts a Matrix and a scalar (scalar on the left)
        /// @param[in] scalar double value to subtract with
        /// @return subtraction of Matrix and scalar
        friend Matrix operator-(double scalar, const Matrix& m);

        /// @brief Transposes matrix
        /// @return Transposed matrix
        Matrix transpose() const;

        /// @brief Inverts matrix
        /// @return Inverted matrix
        Matrix inv();

        /// @brief Extracts subvector from Matrix
        /// @param[in] from starting subvector index
        /// @param[in] to ending subvector index
        /// @return Matrix representing subvector in Matrix
        Matrix subvector(int from, int to) const;

        /// @brief Concatenats Matrixes horizontally
        /// @param[in] other Matrix to concat
        /// @return concatenation of both Matrixes
        Matrix concatRow(const Matrix& other) const;

        /// @brief Extracts column from Matrix
        /// @param[in] column column index
        /// @return Matrix representing a column of the original
        Matrix getColumn(int column) const;

        /// @brief Creates an identity matrix of given size
        /// @param[in] n Size of the identity matrix (n x n)
        /// @return Identity matrix of size n x n
        static Matrix eye(int n);

        /// @brief Calculates Euclidean norm for Matrix
        /// @param[in] m Matrix to calculate Euclidean norm from
        /// @return Returns Euclidean norm for Matrix
        static double norm(const Matrix& m);

        /// @brief Calculates Scalar product of two Matrix
        /// @param[in] m1 Matrix to calculate Scalar product from
        /// @param[in] m2 Matrix to calculate Scalar product from
        /// @return Returns Scalar product of two Matrix
        static double dot(const Matrix& m1, const Matrix& m2);

        /// @brief Calculates Cross product of two Matrix
        /// @param[in] m1 Matrix to calculate Cross product from
        /// @param[in] m2 Matrix to calculate Cross product from
        /// @return Returns Cross product of two Matrix
        static Matrix cross(const Matrix& m1, const Matrix& m2);

        /// @brief Creates a matrix representing a range of integers
        /// @param[in] from starting value of the range
        /// @param[in] to upper bound of the range (exclusive)
        /// @param[in] step step size between values in the range
        /// @return Matrix containing generated sequence
        static Matrix range(int from, int to, int step);

        /// @brief Prints the matrix to the standard output
        void print();

 
    private:
        /// @brief Initializes the internal matrix memory and sets values to zero
        void initMatrix();
 
    private:
        int fil;
        int col;
        double **matrix;
};

#endif
