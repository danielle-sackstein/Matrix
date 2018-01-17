#include <vector>
#include <stdexcept>
#include <iostream>
#include "MatrixException.h"

#pragma once

using std::vector;

template<class T>
class Matrix
{
public:
	Matrix();

	Matrix(unsigned int rows, unsigned int cols);

	Matrix(const Matrix<T> &);

	Matrix(const Matrix<T> &&);

	Matrix(unsigned int rows, unsigned int cols, const vector<T> &cells);

	~Matrix();

	Matrix &operator=(const Matrix<T> &);

	Matrix &operator+(const Matrix<T> &) const;

	Matrix &operator-(const Matrix<T> &) const;

	Matrix &operator*(const Matrix<T> &) const;

	bool operator==(const Matrix<T> &) const;

	bool operator!=(const Matrix<T> &) const;

	bool isSquareMatrix() const;

	Matrix &trans();

	const Matrix &operator()(unsigned int, unsigned int) const;

	Matrix &operator()(unsigned int, unsigned int);

	friend std::ostream &operator<<(std::ostream &, const Matrix<T> &);
	inline unsigned int rows()
	{ return _rows; }

	inline unsigned int cols()
	{ return _cols; }


private:
	unsigned int _hasValidRows(unsigned int rows, unsigned int cols) const;

	unsigned int _hasValidCols(unsigned int rows, unsigned int cols) const;

	unsigned int _rows;
	unsigned int _cols;
	vector<T> _data;

};

template<class T>
Matrix<T>::Matrix()
{
	_rows = 1;
	_cols = 1;
	_data = T{0};
}

template<class T>
Matrix<T>::Matrix(unsigned int rows, unsigned int cols):
		_rows(_hasValidRows(rows, cols)),
		_cols(_hasValidCols(rows, cols)),
		_data(vector<T>(rows * cols, T{0}))
{}

template<class T>
Matrix<T>::Matrix(const Matrix<T> &matrix):
		_rows(matrix._rows),
		_cols(matrix._cols),
		_data(matrix._data) // vector implements operator=
{}

template<class T>
Matrix<T>::Matrix(const Matrix<T> &&matrix):
		_rows(matrix._rows), _cols(matrix._cols), _data(std::move(matrix))
{}

template<class T>
Matrix<T>::Matrix(unsigned int rows, unsigned int cols, const vector<T> &cells):
		_rows(_hasValidRows(rows, cols)),
		_cols(_hasValidCols(rows, cols) ? cols : throw MatrixException("invalid rows or cols")),
		_data(cells)
{}

template<class T>
Matrix<T>::~Matrix()
{
}

template<class T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &matrix)
{
	if (this != &matrix)
	{
		_rows = matrix._rows;
		_cols = matrix._cols;

		//vector implements an operator=
		_data = matrix._data;
	}
	return *this;
}

template<class T>
Matrix<T> &Matrix<T>::operator+(const Matrix<T> &rhs) const
{
	if (_rows != rhs._rows || _cols != rhs._cols)
	{
		throw MatrixException(" ");
	}

	Matrix<T> *matrix = new Matrix<T>(_rows, _cols);

	for (unsigned int i = 0; i < _rows * _cols; i++)
	{
		matrix->_data[i] = matrix->_data[i] + rhs._data[i];
	}
	return *matrix;
}

template<class T>
Matrix<T> &Matrix<T>::operator-(const Matrix<T> &rhs) const
{
	if (_rows != rhs._rows || _cols != rhs._cols)
	{
		throw MatrixException("cannot perform - with different matrices");
	}

	Matrix<T> *matrix = new Matrix<T>(_rows, _cols);
	for (unsigned int i = 0; i < _rows * _cols; i++)
	{
		matrix->_data[i] = matrix->_data[i] - rhs._data[i];
	}
	return *matrix;
}

template<class T>
bool Matrix<T>::operator==(const Matrix<T> &rhs) const
{
	if (_rows != rhs._rows || _cols != rhs._cols)
	{
		throw MatrixException("  ");
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
Matrix<T> &Matrix<T>::operator*(const Matrix<T> &rhs) const
{
	Matrix<T> *matrix = new Matrix<T>(_rows, rhs._cols);
	for (unsigned int i = 0; i < _rows; i++)
	{
		for (unsigned int j = 0; j < rhs._cols; j++)
		{
			unsigned int sum = 0;

			for (unsigned int k = 0; k < rhs._cols; k++)
			{
				sum = sum + this->operator()(i, k) * rhs.operator()(k, j);
			}
			matrix->operator()(i, j) = sum;
		}
	}
	return *matrix;
}

template<class T>
bool Matrix<T>::isSquareMatrix() const
{
	return _rows == _cols;
}

template<class T>
Matrix<T> &Matrix<T>::trans()
{
	if (!this->isSquareMatrix())
	{
		throw MatrixException("cannot perform trans on non squared matrix");
	}

	Matrix<T> *matrix = new Matrix<T>(_rows, _cols);

	for (unsigned int i = 0; i < _rows; i++)
	{
		for (unsigned int j = 0; i < _cols; j++)
		{
			if (i != j)
			{
				matrix->operator()(i, j) = (*this)(j, i);
			}
		}
	}
	return *matrix;
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
Matrix<T> &Matrix<T>::operator()(unsigned int row, unsigned int col)
{
	return _data[row * col + col];
}

template<class T>
const Matrix<T> &Matrix<T>::operator()(unsigned int row, unsigned int col) const
{
	return _data[row * col + col];
}

template<class T>
unsigned int Matrix<T>::_hasValidCols(unsigned int rows, unsigned int cols) const
{
	if ((rows == 0 && cols != 0) || (rows != 0 && cols == 0))
	{
		throw MatrixException("invalid rows or cols");
	}
	return _cols;
}

template<class T>
unsigned int Matrix<T>::_hasValidRows(unsigned int rows, unsigned int cols) const
{
	if ((rows == 0 && cols != 0) || (rows != 0 && cols == 0))
	{
		throw MatrixException("invalid rows or cols");
	}
	return _rows;
}
