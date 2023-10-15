#include <iostream>
#include <random> // generate random number

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
    void print()
    {
        for (int i = 0; i < _rows; i++) {
            for (int j = 0; j < _cols; j++) {
                std::cout << _data[i][j] << " ";
            }
            std::cout << std::endl;
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


    T& operator()(size_t r, size_t c)
    {
        return _data[r][c];
    }

    Matrix operator+(const Matrix& other) const
    {
        Matrix c(_rows, _cols, 0);
        for(size_t i = 0; i < _rows; i++)
        {
            for(size_t j = 0; j < _cols; j++)
            {
                c._data[i][j] = _data[i][j] + other._data[i][j];
            }
        }
        return c;
    }
    Matrix operator-(const Matrix& other);

    T getRows()
    {
        return _rows;
    }

    T getCols()
    {
        return _cols;
    }
};

int main()
{
    Matrix a(2, 2, 4);
    a.print();
    std::cout << std::endl;
    Matrix b(2, 2, 1, 3);
    b.print();
    std::cout << std::endl;
    Matrix c = a + b;
    c.print();
}