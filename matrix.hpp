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

	Matrix &operator=(const Matrix &);

	Matrix &operator+(const Matrix &) const;

	Matrix &operator-(const Matrix &) const;

	Matrix &operator*(const Matrix &) const;

	bool operator==(const Matrix &) const;

	bool operator!=(const Matrix &) const;

	bool isSquareMatrix() const;

	Matrix &trans();

	const Matrix &operator()(unsigned int, unsigned int) const;

	Matrix &operator()(unsigned int, unsigned int);

	friend std::ostream &operator<<(std::ostream &, const Matrix &);

	inline unsigned int rows()
	{ return _rows; }

	inline unsigned int cols()
	{ return _cols; }


private:
	bool _hasValidCords(unsigned int rows, unsigned int cols) const;

	unsigned int _rows;
	unsigned int _cols;
	vector _data;

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
		_rows(_hasValidCords(rows, cols) ? rows : throw new MatrixException),
		_cols(_hasValidCords(rows, cols) ? cols : throw new MatrixException),
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
		_rows(_hasValidCords(rows, cols) ? rows : throw new MatrixException),
		_cols(_hasValidCords(rows, cols) ? cols : throw new MatrixException),
		_data(cells)
{}

template<class T>
Matrix<T>::~Matrix()
{
	delete Matrix<T>;
}

template<class T>
Matrix &Matrix<T>::operator=(const Matrix &matrix)
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
Matrix &Matrix<T>::operator+(const Matrix &rhs) const
{
	if (_rows != rhs._rows || _cols != rhs._cols)
	{
		throw new MatrixException;
	}

	Matrix<T> *matrix = new Matrix<T>(_rows, _cols);

	for (unsigned int i = 0; i < _rows * _cols; i++)
	{
		matrix->_data[i] = matrix->_data[i] + rhs._data[i];
	}
	return *matrix;
}

template<class T>
Matrix &Matrix<T>::operator-(const Matrix &rhs) const
{
	if (_rows != rhs._rows || _cols != rhs._cols)
	{
		throw new MatrixException;
	}

	Matrix<T> *matrix = new Matrix<T>(_rows, _cols);
	for (unsigned int i = 0; i < _rows * _cols; i++)
	{
		matrix->_data[i] = matrix->_data[i] - rhs._data[i];
	}
	return *matrix;
}

template<class T>
bool Matrix<T>::operator==(const Matrix &rhs) const
{
	if (_rows != rhs._rows || _cols != rhs._cols)
	{
		throw new MatrixException;
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
bool Matrix<T>::operator!=(const Matrix &rhs) const
{
	return !operator==(rhs);
}

template<class T>
Matrix &Matrix::operator*(const Matrix &rhs) const
{
	Matrix<T> *matrix = new Matrix<T>(_rows, rhs._cols);
	for (unsigned int i = 0; i < _rows; i++)
	{
		for (unsigned int j = 0; j < rhs._cols; j++)
		{
			unsigned int sum = 0;

			for (unsigned int k = 0; k < rhs._cols; k++)
			{
				sum = sum + _data(i, k) * rhs._data(k, j);
			}
			matrix->(i, j) = sum;
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
Matrix &Matrix<T>::trans()
{
	if (!this->isSquareMatrix())
	{
		throw new MatrixException;
	}

	Matrix<T> *matrix = new Matrix<T>(_rows, _cols);

	for (unsigned int i = 0; i < _rows; i++)
	{
		for (unsigned int j = 0; i < _cols; j++)
		{
			if (i != j)
			{
				matrix->(i, j) = this->(j, i);
			}
		}
	}
	return *matrix;
}


template<class T>
friend std::ostream &Matrix<T>::operator<<(std::ostream &os, const Matrix &matrix)
{
	for (unsigned int i = 0; i < matrix._rows; i++)
	{
		for (unsigned int j = 0; j < matrix._rows; j++)
		{
			os << matrix(i, j) << "	";
		}
		os << std::endl;
	}
	return os;
}

template<class T>
Matrix &Matrix<T>::operator()(unsigned int row, unsigned int col)
{
	return _data[row * col + col];
}

template<class T>
const Matrix &Matrix<T>::operator()(unsigned int row, unsigned int col) const
{
	return _data[row * col + col];
}

bool _hasValidCords(unsigned int rows, unsigned int cols) const
{
	return !(rows == 0 && cols != 0) && !(rows != 0 && cols == 0);
}
