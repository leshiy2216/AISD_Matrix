#include <iostream>
#include <random> // generate random number
#include <complex>
#include <stdexcept>

template<class T>
class Matrix
{
private:
    size_t _rows;
    size_t _cols;
    T** _data;
public:
    Matrix(size_t r, size_t c, T value)
    {
        this->_rows = r;
        this->_cols = c;
        _data = new T * [_rows];
        for (size_t i = 0; i < _rows; i++)
        {
            _data[i] = new T[_cols];
            for (size_t j = 0; j < _cols; j++)
            {
                _data[i][j] = value;
            }
        }
    }

    Matrix(size_t r, size_t c, T low, T up)
    {
        this->_rows = r;
        this->_cols = c;
        _data = new T * [_rows];
        std::mt19937 rng(std::random_device{}());
        std::uniform_real_distribution<T> distribution(low, up);
        for (size_t i = 0; i < _rows; i++)
        {
            _data[i] = new T[_cols];
            for (size_t j = 0; j < _cols; j++)
            {
                _data[i][j] = distribution(rng);
            }
        }
    }

    ~Matrix()
    {
        for (size_t i = 0; i < _rows; i++)
        {
            delete[] _data[i];
        }
        delete[] _data;
    }

    void print()
    {
        for (size_t i = 0; i < _rows; i++)
        {
            for (size_t j = 0; j < _cols; j++)
            {
                std::cout << _data[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    T& operator()(size_t r, size_t c) const
    {
        return _data[r][c];
    }

    Matrix operator+(const Matrix& other) const
    {
        Matrix c(_rows, _cols, 0);
        if ((_rows != other.getRows()) || (_cols != other.getCols())) throw std::runtime_error("Matrix's can be equally");

        for (size_t i = 0; i < _rows; i++)
        {
            for (size_t j = 0; j < _cols; j++)
            {
                c(i, j) = _data[i][j] + other(i, j);
            }
        }
        return c;
    }

    Matrix operator-(const Matrix& other) const
    {
        Matrix c(_rows, _cols, 0);
        if ((_rows != other.getRows()) || (_cols != other.getCols())) throw std::runtime_error("Matrix's can be equally");

        for (size_t i = 0; i < _rows; i++)
        {
            for (size_t j = 0; j < _cols; j++)
            {
                c(i, j) = _data[i][j] - other(i, j);
            }
        }
        return c;
    }

    Matrix operator*(const Matrix& other) const
    {
        if (_cols != other.getRows()) throw std::runtime_error("Multiply is impossible");
        Matrix c(_rows, other.getCols(), 0);
        for (size_t i = 0; i < _rows; i++)
        {
            for (size_t j = 0; j < other.getCols(); j++)
            {
                for (size_t k = 0; k < _cols; k++)
                {
                    c(i, j) += _data[i][k] * other(k, j);
                }
            }
        }
        return c;
    }

    Matrix operator*(T scalar) const
    {
        Matrix c(_rows, _cols, 0);
        for (size_t i = 0; i < _rows; i++)
        {
            for (size_t j = 0; j < _cols; j++)
            {
                c(i, j) = _data[i][j] * scalar;
            }
        }
        return c;
    }

    friend Matrix operator*(T scalar, const Matrix& other)
    {
        Matrix c(other.getRows(), other.getCols(), 0);
        for (size_t i = 0; i < other.getRows(); i++)
        {
            for (size_t j = 0; j < other.getCols(); j++)
            {
                c(i, j) = other(i, j) * scalar;
            }
        }
        return c;
    }

    Matrix operator/(T scalar) const
    {
        if (scalar == 0) throw std::runtime_error("Division by zero");

        Matrix c(_rows, _cols, 0);
        for (size_t i = 0; i < _rows; ++i)
        {
            for (size_t j = 0; j < _cols; ++j)
            {
                c(i, j) = _data[i][j] / scalar;
            }
        }
        return c;
    }

    T trace() const
    {
        T sum = 0;
        if (_rows != _cols) throw std::runtime_error("Only for square matrix's");

        for (size_t i = 0; i < _rows; i++)
        {
            sum = sum + _data[i][i];
        }
        return sum;
    }

    T determinant() const
    {
        if (_rows != _cols) throw std::runtime_error("Determinant only for square matrix's");
        if (_rows == 1) return _data[0][0];
        if (_rows == 2) return _data[0][0] * _data[1][1] - _data[0][1] * _data[1][0];
        else
        {
            T det = 0;
            for (size_t j = 0; j < _cols; j++)
            {
                Matrix submatrix(_rows - 1, _cols - 1, 0);
                for (size_t row = 1; row < _rows; row++)
                {
                    for (size_t col = 0; col < _cols; col++)
                    {
                        if (col == j)
                            continue;
                        submatrix(row - 1, col < j ? col : col - 1) = _data[row][col];
                    }
                }
                det += (_data[0][j] * submatrix.determinant()) * (j % 2 == 0 ? 1 : -1);
            }
            return det;
        }
    }

    Matrix reverseMatrix()
    {
        if (_rows != _cols) throw std::runtime_error("Reverse only for square matrix's");
        T det = determinant();
        if (det == 0) throw std::runtime_error("Cannot calculate reverse matrix");

        Matrix reversed(_rows, _rows, 0);
        for (size_t i = 0; i < _rows; i++)
        {
            for (size_t j = 0; j < _rows; j++)
            {
                Matrix minor(_rows - 1, _rows - 1, 0);
                for (size_t row = 0, minorRow = 0; row < _rows; row++)
                {
                    if (row == i) continue;
                    for (size_t col = 0, minorCol = 0; col < _rows; col++)
                    {
                        if (col == j) continue;
                        minor(minorRow, minorCol) = _data[row][col];
                        minorCol++;
                    }
                    minorRow++;
                }
                T cofactor = minor.determinant() * ((i + j) % 2 == 0 ? 1 : -1);
                reversed(j, i) = cofactor / det;
            }
        }
        return reversed;
    }

    size_t getRows() const
    {
        return _rows;
    }

    size_t getCols() const
    {
        return _cols;
    }
};

template <class T>
class Matrix<std::complex<T>>
{
private:
    size_t _rows;
    size_t _cols;
    std::complex<T>** _data;

public:
    Matrix(size_t r, size_t c, std::complex<T> value)
        : _rows(r), _cols(c)
    {
        _data = new std::complex<T>*[_rows];
        for (size_t i = 0; i < _rows; i++)
        {
            _data[i] = new std::complex<T>[_cols];
            for (size_t j = 0; j < _cols; j++)
            {
                _data[i][j] = value;
            }
        }
    }

    Matrix(size_t r, size_t c, std::complex<T> low, std::complex<T> up)
        : _rows(r), _cols(c)
    {
        _data = new std::complex<T>*[_rows];
        std::mt19937 rng(std::random_device{}());
        std::uniform_real_distribution<T> realDistribution(low.real(), up.real());
        std::uniform_real_distribution<T> imagDistribution(low.imag(), up.imag());
        for (size_t i = 0; i < _rows; i++)
        {
            _data[i] = new std::complex<T>[_cols];
            for (size_t j = 0; j < _cols; j++)
            {
                _data[i][j] = std::complex<T>(realDistribution(rng), imagDistribution(rng));
            }
        }
    }

    ~Matrix()
    {
        for (size_t i = 0; i < _rows; i++)
        {
            delete[] _data[i];
        }
        delete[] _data;
    }

    void print()
    {
        for (size_t i = 0; i < _rows; i++)
        {
            for (size_t j = 0; j < _cols; j++)
            {
                std::cout << _data[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    Matrix<std::complex<T>> operator+(const Matrix<std::complex<T>>& other) const
    {
        Matrix<std::complex<T>> result(_rows, _cols, std::complex<T>(0, 0));
        if ((_rows != other.getRows()) || (_cols != other.getCols())) throw std::runtime_error("Матрицы должны быть одинакового размера");

        for (size_t i = 0; i < _rows; i++)
        {
            for (size_t j = 0; j < _cols; j++)
            {
                result(i, j) = _data[i][j] + other(i, j);
            }
        }
        return result;
    }

    Matrix<std::complex<T>> operator-(const Matrix<std::complex<T>>& other) const
    {
        Matrix<std::complex<T>> result(_rows, _cols, std::complex<T>(0, 0));
        if ((_rows != other.getRows()) || (_cols != other.getCols())) throw std::runtime_error("Матрицы должны быть одинакового размера");

        for (size_t i = 0; i < _rows; i++)
        {
            for (size_t j = 0; j < _cols; j++)
            {
                result(i, j) = _data[i][j] - other(i, j);
            }
        }
        return result;
    }

    Matrix<std::complex<T>> operator*(const Matrix<std::complex<T>>& other) const
    {
        if (_cols != other.getRows()) throw std::runtime_error("Умножение невозможно");
        Matrix<std::complex<T>> result(_rows, other.getCols(), std::complex<T>(0, 0));

        for (size_t i = 0; i < _rows; i++)
        {
            for (size_t j = 0; j < other.getCols(); j++)
            {
                for (size_t k = 0; k < _cols; k++)
                {
                    result(i, j) += _data[i][k] * other(k, j);
                }
            }
        }
        return result;
    }

    Matrix<std::complex<T>> operator*(std::complex<T> scalar) const
    {
        Matrix<std::complex<T>> result(_rows, _cols, std::complex<T>(0, 0));
        for (size_t i = 0; i < _rows; i++)
        {
            for (size_t j = 0; j < _cols; j++)
            {
                result(i, j) = _data[i][j] * scalar;
            }
        }
        return result;
    }

    friend Matrix<std::complex<T>> operator*(std::complex<T> scalar, const Matrix<std::complex<T>>& matrix)
    {
        Matrix<std::complex<T>> result(matrix.getRows(), matrix.getCols(), std::complex<T>(0, 0));
        for (size_t i = 0; i < matrix.getRows(); i++)
        {
            for (size_t j = 0; j < matrix.getCols(); j++)
            {
                result(i, j) = matrix(i, j) * scalar;
            }
        }
        return result;
    }

    Matrix<std::complex<T>> operator/(std::complex<T> scalar) const
    {
        if (scalar == std::complex<T>(0, 0)) throw std::runtime_error("Деление на ноль невозможно");

        Matrix<std::complex<T>> result(_rows, _cols, std::complex<T>(0, 0));
        for (size_t i = 0; i < _rows; ++i)
        {
            for (size_t j = 0; j < _cols; ++j)
            {
                result(i, j) = _data[i][j] / scalar;
            }
        }
        return result;
    }

    std::complex<T> trace() const
    {
        std::complex<T> sum(0, 0);
        if (_rows != _cols) throw std::runtime_error("Только для квадратных матриц");

        for (size_t i = 0; i < _rows; i++)
        {
            sum = sum + _data[i][i];
        }
        return sum;
    }

    std::complex<T> determinant() const
    {
        if (_rows != _cols) throw std::runtime_error("Определитель только для квадратных матриц");
        if (_rows == 1) return _data[0][0];
        if (_rows == 2) return _data[0][0] * _data[1][1] - _data[0][1] * _data[1][0];
        else
        {
            std::complex<T> det = std::complex<T>(0, 0);
            for (size_t j = 0; j < _cols; j++)
            {
                Matrix<std::complex<T>> submatrix(_rows - 1, _cols - 1, std::complex<T>(0, 0));
                for (size_t row = 1; row < _rows; row++)
                {
                    for (size_t col = 0; col < _cols; col++)
                    {
                        if (col == j)
                            continue;
                        submatrix(row - 1, col < j ? col : col - 1) = _data[row][col];
                    }
                }
                det += (_data[0][j] * submatrix.determinant()) * (j % 2 == 0 ? std::complex<T>(1, 0) : std::complex<T>(-1, 0));
            }
            return det;
        }
    }

    Matrix<std::complex<T>> reverseMatrix()
    {
        if (_rows != _cols) throw std::runtime_error("Обратная матрица только для квадратных матриц");
        std::complex<T> det = determinant();
        if (det == std::complex<T>(0, 0)) throw std::runtime_error("Невозможно вычислить обратную матрицу");

        Matrix<std::complex<T>> reversed(_rows, _rows, std::complex<T>(0, 0));
        for (size_t i = 0; i < _rows; i++)
        {
            for (size_t j = 0; j < _rows; j++)
            {
                Matrix<std::complex<T>> minor(_rows - 1, _rows - 1, std::complex<T>(0, 0));
                for (size_t row = 0, minorRow = 0; row < _rows; row++)
                {
                    if (row == i) continue;
                    for (size_t col = 0, minorCol = 0; col < _rows; col++)
                    {
                        if (col == j) continue;
                        minor(minorRow, minorCol) = _data[row][col];
                        minorCol++;
                    }
                    minorRow++;
                }
                std::complex<T> cofactor = minor.determinant() * ((i + j) % 2 == 0 ? std::complex<T>(1, 0) : std::complex<T>(-1, 0));
                reversed(j, i) = cofactor / det;
            }
        }
        return reversed;
    }

    size_t getRows() const
    {
        return _rows;
    }

    size_t getCols() const
    {
        return _cols;
    }

    std::complex<T>& operator()(size_t r, size_t c) const
    {
        return _data[r][c];
    }
};

int main()
{
    setlocale(LC_ALL, "Rus");

    Matrix<double> A(3, 3, 1, 5);
    Matrix<double> b(3, 1, 1, 5);
    b(0, 0) = 2;
    b(1, 0) = 4;
    b(2, 0) = 6;

    std::cout << "Матрица A:\n";
    A.print();
    std::cout << "\n";

    std::cout << "Вектор b:\n";
    b.print();
    std::cout << "\n";

    Matrix<double> x = A.reverseMatrix() * b;

    std::cout << "Решение уравнения Ax = b:\n";
    x.print();
    std::cout << "\n";

    return 0;
}