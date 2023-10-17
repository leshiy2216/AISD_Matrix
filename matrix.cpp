#include <iostream>
#include <random> // generate random number
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
        _data = new T*[_rows];
        for(size_t i = 0; i < _rows; i++)
        {
            _data[i] = new T[_cols];
            for(size_t j = 0; j < _cols; j++)
            {
                _data[i][j] = value;
            }
        }
    }
    
    Matrix(size_t r, size_t c, T low, T up)
    {
        this->_rows = r;
        this->_cols = c;
        _data = new T*[_rows];
        std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<T> distribution(low, up);
        for(size_t i = 0; i < _rows; i++)
        {
            _data[i] = new T[_cols];
            for(size_t j = 0; j < _cols; j++)
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

        for(size_t i = 0; i < _rows; i++)
        {
            for(size_t j = 0; j < _cols; j++)
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

        for(size_t i = 0; i < _rows; i++)
        {
            for(size_t j = 0; j < _cols; j++)
            {
                c(i, j) = _data[i][j] - other(i, j);
            }
        }
        return c;
    }

    Matrix operator*(const Matrix& other) const
    {
        if (_cols!=other.getRows()) throw std::runtime_error("Multiply is impossible");
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
        for(size_t i = 0; i < _rows; i++)
        {
            for(size_t j = 0; j < _cols; j++)
            {
                c(i, j) = _data[i][j] * scalar;
            }
        }
        return c;
    }

    friend Matrix operator*(T scalar, const Matrix& other)
    {
        Matrix c(other.getRows(), other.getCols(), 0);
        for(size_t i = 0; i < other.getRows(); i++)
        {
            for(size_t j = 0; j < other.getCols(); j++)
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

int main()
{
    Matrix a(3, 3, 1, 3);
    a.print();
    std::cout << std::endl;
    Matrix b(3, 3, 3, 8);
    b.print();
    Matrix x = a.reverseMatrix() * b;
    std::cout << std::endl;
    x.print();
}