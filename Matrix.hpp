#include <vector>
#include <exception>
#include <iostream>
#include <thread>

#include "MatrixException.h"
#include "Complex.h"

#pragma once

template<class T>
class Matrix
{
public:

	Matrix();

	Matrix(unsigned int rows, unsigned int cols);

	Matrix(const Matrix<T> &);

	Matrix(Matrix<T> &&);

	Matrix(unsigned int rows, unsigned int cols, const std::vector<T> &cells);

	~Matrix();

	Matrix<T> &operator=(const Matrix<T> &);

	Matrix<T> operator+(const Matrix<T> &) const;

	Matrix<T> operator-(const Matrix<T> &) const;

	Matrix<T> operator*(const Matrix<T> &) const;

	bool operator==(const Matrix<T> &) const;

	bool operator!=(const Matrix<T> &) const;

	bool isSquareMatrix() const;

	Matrix<T> trans() const;

	const T &operator()(unsigned int, unsigned int) const;

	T &operator()(unsigned int, unsigned int);

	template<typename U>
	friend std::ostream &operator<<(std::ostream &, const Matrix<U> &);

	inline unsigned int rows() const
	{ return _rows; }

	inline unsigned int cols() const
	{ return _cols; }

	typedef typename std::vector<T>::const_iterator const_iterator;

	const_iterator begin() const
	{ return _data.begin(); }

	const_iterator end() const
	{ return _data.begin(); }

    // TODO: Print
    static void setParallel(bool isParallel) { _isParallel = isParallel; }

private:

    static Matrix<T> multSerial(const Matrix<T> &lhs, const Matrix<T> &rhs);
    static Matrix<T> multParallel(const Matrix<T> &lhs, const Matrix<T> &rhs);

	unsigned int _rows;
	unsigned int _cols;
	std::vector<T> _data;
    static bool _isParallel;
};

template<typename T>
bool Matrix<T>::_isParallel = true;

template<class T>
Matrix<T>::Matrix() :
		_rows(1),
		_cols(1),
		_data(1, 0)
{
}

template<class T>
Matrix<T>::Matrix(unsigned int rows, unsigned int cols):
		Matrix(rows, cols, std::vector<T>(rows * cols))
{
}

template<class T>
Matrix<T>::Matrix(const Matrix<T> &matrix):
		_rows(matrix._rows),
		_cols(matrix._cols),
		_data(matrix._data)
{}

template<class T>
Matrix<T>::Matrix(Matrix<T> &&matrix):
		_rows(matrix._rows),
		_cols(matrix._cols),
		_data(std::move(matrix._data))
{}

template<class T>
Matrix<T>::Matrix(unsigned int rows, unsigned int cols, const std::vector<T> &cells):
		_rows(rows),
		_cols(cols),
		_data(cells)
{
	if ((rows == 0) || (cols == 0))
	{
		throw MatrixException("invalid dimensions\n");
	}
}

template<class T>
Matrix<T>::~Matrix()
{}

template<class T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &matrix)
{
	if (this != &matrix)
	{
		_data = matrix._data;

		_rows = matrix._rows;
		_cols = matrix._cols;
	}
	return *this;
}

template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &rhs) const
{
	if ((_rows != rhs._rows) || (_cols != rhs._cols))
	{
		throw MatrixException("Incompatible dimensions");
	}

	Matrix<T> matrix(_rows, _cols);

	for (unsigned int i = 0; i < _rows * _cols; i++)
	{
		matrix._data[i] = _data[i] + rhs._data[i];
	}

	return matrix;
}

template<class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &rhs) const
{
	if (_rows != rhs._rows || _cols != rhs._cols)
	{
		throw MatrixException("Incompatible dimensions");
	}

	Matrix<T> matrix(_rows, _cols);

	for (unsigned int i = 0; i < _rows * _cols; i++)
	{
		matrix._data[i] = _data[i] - rhs._data[i];
	}

	return matrix;
}

template<class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &rhs) const
{
	if (_cols != rhs._rows)
	{
		throw MatrixException("Incompatible dimensions");
	}

	return _isParallel
		   ? multParallel(*this, rhs)
		   : multSerial(*this, rhs);
}

template<class T>
Matrix<T> Matrix<T>::multSerial(const Matrix<T> &lhs, const Matrix<T> &rhs)
{
    Matrix<T> result (lhs._rows, rhs._cols);

    for (unsigned int i = 0; i < lhs._rows; i++)
    {
        for (unsigned int j = 0; j < rhs._cols; j++)
        {
            T sum = T{};

            for (unsigned int k = 0; k < rhs._cols; k++)
            {
                sum = sum + lhs(i, k) * rhs(k, j);
            }

            result(i, j) = sum;
        }
    }

    return result;
}

template <typename T>
class VectorMultiplier
{
public:

    VectorMultiplier(const Matrix<T> &lhs, const Matrix<T> &rhs, Matrix<T>& result) :
        _lhs(lhs),
        _rhs(rhs),
        _result(result)
    {}

    void operator()(unsigned int i)
    {
        for (unsigned int j = 0; j < _rhs.cols(); j++)
        {
            T sum = T{};

            for (unsigned int k = 0; k < _rhs.cols(); k++)
            {
                sum = sum + _lhs(i, k) * _rhs(k, j);
            }

            _result(i, j) = sum;
        }
    }

private:

    const Matrix<T> &_lhs;
    const Matrix<T> &_rhs;
    Matrix<T>& _result;
};

template<class T>
Matrix<T> Matrix<T>::multParallel(const Matrix<T> &lhs, const Matrix<T> &rhs)
{
    Matrix<T> result (lhs._rows, rhs._cols);

    VectorMultiplier<T> multiplier (lhs, rhs, result);

    std::vector<std::thread> threads (lhs._rows);

    for (unsigned int i = 0; i < lhs._rows; i++)
    {
        threads[i] = std::thread (std::ref(multiplier), i);
    }

    for (auto& t : threads){
        t.join();
    }

    return result;
}

template<class T>
bool Matrix<T>::operator==(const Matrix<T> &rhs) const
{
	if ((_rows != rhs._rows) || (_cols != rhs._cols))
	{
		return false;
	}

	for (unsigned int i = 0; i < _rows * _cols; i++)
	{
		if (_data[i] != rhs._data[i])
		{
			return false;
		}
	}
	return true;
}

template<class T>
bool Matrix<T>::operator!=(const Matrix<T> &rhs) const
{
	return !operator==(rhs);
}

template<class T>
bool Matrix<T>::isSquareMatrix() const
{
	return _rows == _cols;
}

template<class T>
Matrix<T> Matrix<T>::trans() const
{
	if (!isSquareMatrix())
	{
		throw MatrixException("Cannot perform trans on non squared matrix\n");
	}

	Matrix<T> matrix(_rows, _cols);
	const Matrix<T> &_this = *this;

	for (unsigned int i = 0; i < _rows; i++)
	{
		for (unsigned int j = 0; j < _cols; j++)
		{
			matrix(i, j) = _this(j, i);
		}
	}
	return matrix;
}

template<>
inline Matrix<Complex> Matrix<Complex>::trans() const
{
	if (!isSquareMatrix())
	{
		throw MatrixException("Cannot perform trans on non squared matrix\n");
	}

	Matrix<Complex> matrix(_rows, _cols);
	const Matrix<Complex> &_this = *this;

	for (unsigned int i = 0; i < _rows; ++i)
	{
		for (unsigned int j = 0; j < _cols; ++j)
		{
			matrix(j, i) = _this(i, j).conj();
		}
	}
	return matrix;
}

template<class T>
std::ostream &operator<<(std::ostream &os, const Matrix<T> &matrix)
{
	for (unsigned int i = 0; i < matrix._rows; i++)
	{
		for (unsigned int j = 0; j < matrix._cols; j++)
		{
			os << matrix(i, j) << "	";
		}
		os << std::endl;
	}
	return os;
}

template<class T>
T &Matrix<T>::operator()(unsigned int row, unsigned int col)
{
	if ((row > _rows) || (row < 0))
	{
		throw MatrixException("Row out of bounds");
	}
	if ((col > _cols) || (col < 0))
	{
		throw MatrixException("Col out of bounds");
	}
	return _data[row * _cols + col];
}

template<class T>
const T &Matrix<T>::operator()(unsigned int row, unsigned int col) const
{
	if ((row > _rows) || (row < 0))
	{
		throw MatrixException("Row out of bounds");
	}
	if ((col > _cols) || (col < 0))
	{
		throw MatrixException("Col out of bounds");
	}
	return _data[row * _cols + col];
}


