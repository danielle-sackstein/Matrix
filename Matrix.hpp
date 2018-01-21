#include <vector>
#include <exception>
#include <iostream>
#include <thread>

#include "MatrixException.h"
#include "Complex.h"

#pragma once

#define ERROR_IN_CREATING_MATRIX "Cannot create a matrix with 0 rows or columns"

#define ERROR_IN_PERFORMING_ADD_OPERATOR "Incompatible dimensions"

#define ERROR_IN_PERFORMING_REDUCTION_OPERATOR "Incompatible dimensions"

#define ERROR_IN_PERFORMING_MULT_OPERATOR "Incompatible dimensions"

#define ERROR_IN_PERFORMING_TRANS "Cannot perform trans on non squared matrix\n"

#define INDEX_OUT_OF_RANGE_OF_ROWS "Row out of bounds"

#define INDEX_OUT_OF_RANGE_OF_COLS "Col out of bounds"

#define PARALLEL_MODE "Parallel"

#define NON_PARALLEL_MODE "non-Parallel"

#define CHANGE_MODE "Generic Matrix mode changed to "

#define CHANGE_MODE_PART_2 " mode"

/**This class is a template library of Matrix. It has various public methods which perform
 * operations on matrices.**/

template<class T>
class Matrix
{
public:

	typedef typename std::vector<T>::const_iterator const_iterator;

	/**
	 * Default constructor
	 */
	Matrix();

	/**
	* A constructor that receives dimensions of a matrix and initializes the matrix with
	* this size while initializing the cells to the zero elements.
	* @param rows - number of rows.
	* @param cols - number of columns.
	*/
	Matrix(unsigned int rows, unsigned int cols);

	/**
   * A copy constructor.
   * @param matrix - the matrix we want to copy.
   */
	Matrix(const Matrix<T> &);

	/**
	* move constructor.
	* @param matrix - the matrix we would like to move.
	*/
	Matrix(Matrix<T> &&);

	/**
     * A constructor that initializes a new matrix with a specific vector.
     * @param rows - number of rows.
     * @param cols - number of columns.
     * @param cells - the vector with the elements in the matrix.
     */
	Matrix(unsigned int rows, unsigned int cols, const std::vector<T> &cells);

	/**
   * A default destructor.
   */
	~Matrix();

	/**
   * Assignment operator
   * @param matrix - the receives matrix.
   * @return this after it was assigned by the receives matrix.
   */
	Matrix<T> &operator=(const Matrix<T> &);

	/**
     * + operator. Receives a matrix and returns a new matrix that represents the addition of
     * the original matrix and the received matrix.
     * @param matrix - the receives matrix.
     * @return a new matrix that is the addition of the original and the receives matrix.
     */
	Matrix<T> operator+(const Matrix<T> &) const;

	/**
     * - operator. Receives a matrix and returns a new matrix that represents the reduction of
     * the original matrix and the received matrix.
     * @param matrix - the receives matrix.
     * @return a new matrix that is the reduction of the original and the receives matrix.
     */
	Matrix<T> operator-(const Matrix<T> &) const;

	/**
	* * operator. Receives a matrix and returns a new matrix that represents the multiplication of
	* the original matrix and the receives matrix that the operator has received.
	* @param matrix - the receives matrix.
	* @return a new matrix that is the product of the original matrix and the receives matrix.
	*/
	Matrix<T> operator*(const Matrix<T> &) const;

	/**
	* This operator receives a matrix and checks whether it is the same as the received matrix.
	* @param matrix - the matrix we are comparing to.
	* @return true if the matrix is equals to this, false otherwise.
	*/
	bool operator==(const Matrix<T> &) const;

	/**
	* This operator receives a matrix and checks whether it is not the same as the received matrix.
	* @param matrix - the matrix we are comparing to.
	* @return true if the matrix is not the same as - this, false otherwise.
	*/
	bool operator!=(const Matrix<T> &) const;

	/**
	* @return true iff this is quadratic, false otherwise.
	*/
	bool isSquareMatrix() const;

	/**
	* trans() operator . Returns a new matrix that is the trans of - this.
	* @return a new matrix that is the trans of - this.
	*/
	Matrix<T> trans() const;

	/**
	* This operator receives a row and a column and returns the const element located in the [row,
	* col] int this matrix.
	* @param row - The row of the element.
	* @param col - The column of the element.
	* @return the element which is located in the given row and column
	*/
	const T &operator()(unsigned int, unsigned int) const;

	/**
	* This operator receives a row and a column and returns the element located in the [row,
	* col] int this matrix.
	* @param row - The row of the element.
	* @param col - The column of the element.
	* @return the element which is located in the given row and column
	*/
	T &operator()(unsigned int, unsigned int);

	/**
     * This friend operator receives a matrix and an output stream and prints the matrix.
     * @tparam P - a class type.
     * @param os - output stream.
     * @param print - matrix to print.
     * @return the received output stream.
     */
	template<typename U>
	friend std::ostream &operator<<(std::ostream &, const Matrix<U> &);

	/**
	 * @return the number of rows
	 */
	inline unsigned int rows() const
	{ return _rows; }

	/**
	 * @return the number of cols
	 */
	inline unsigned int cols() const
	{ return _cols; }

	/**
	* @return a const iterator which can iterate on all matrix elements.
	*/
	const_iterator begin() const
	{ return _data.begin(); }

	/**
     * @return a const iterator that to the end of the matrix(after the last argument).
     */
	const_iterator end() const
	{ return _data.begin(); }

	/**
     * This method receives a bool value which determined the performance of the ", * operations.
     * If is is true, the multiplication and plus operations of all the matrices will be performed
     * by parallel programming.
     * Otherwise (false value) will return the class to its default behavior (withous using
     * parallel programming).
     * @param parallel - a bool value which determines the mode of the operations +, *
     * true if parallel, false if non-parallel.
     */
	static void setParallel(bool isParallel);

private:

	/**
	 * performs a regular multification.
	 * @param lhs left matrix
	 * @param rhs right matrix
	 * @return the new matrix after the progress
	 */
	static Matrix<T> _multSerial(const Matrix<T> &lhs, const Matrix<T> &rhs);

	/**
	 * performs a regular add.
	 * @param lhs left matrix
	 * @param rhs right matrix
	 * @return the new matrix after the progress
	 */
	static Matrix<T> _addSerial(const Matrix<T> &lhs, const Matrix<T> &rhs);

	/**
	 * Performs a parallel programming multification
	 * @param lhs left matrix
	 * @param rhs right matrix
	 * @return the new matrix after the progress
	 */
	static Matrix<T> _multParallel(const Matrix<T> &lhs, const Matrix<T> &rhs);

	/**
	 * Performs a parallel programming add
	 * @param lhs left matrix
	 * @param rhs right matrix
	 * @return the new matrix after the progress
	 */
	static Matrix<T> _addParallel(const Matrix<T> &lhs, const Matrix<T> &rhs);

	unsigned int _rows;
	unsigned int _cols;
	std::vector<T> _data;
	static bool _isParallel;
};

template<typename T>
bool Matrix<T>::_isParallel = true;

/**
* Default constructor
*/
template<class T>
Matrix<T>::Matrix() :
		_rows(1),
		_cols(1),
		_data(1, 0)
{}

/**
* A constructor that receives dimensions of a matrix and initializes the matrix with
* this size while initializing the cells to the zero elements.
* @param rows - number of rows.
* @param cols - number of columns.
*/
template<class T>
Matrix<T>::Matrix(unsigned int rows, unsigned int cols):
		Matrix(rows, cols, std::vector<T>(rows * cols))
{}

/**
* A copy constructor.
* @param matrix - the matrix we want to copy.
*/
template<class T>
Matrix<T>::Matrix(const Matrix<T> &matrix):
		_rows(matrix._rows),
		_cols(matrix._cols),
		_data(matrix._data)
{}

/**
* move constructor.
* @param matrix - the matrix we would like to move.
*/
template<class T>
Matrix<T>::Matrix(Matrix<T> &&matrix):
		_rows(matrix._rows),
		_cols(matrix._cols),
		_data(std::move(matrix._data))
{}

/**
 * A constructor that initializes a new matrix with a specific vector.
 * @param rows - number of rows.
 * @param cols - number of columns.
 * @param cells - the vector with the elements in the matrix.
 */
//TODO check about exception
template<class T>
Matrix<T>::Matrix(unsigned int rows, unsigned int cols, const std::vector<T> &cells):
		_rows(rows),
		_cols(cols),
		_data(cells)
{
	if ((rows == 0) || (cols == 0) || cells.size() != _rows*_cols)
	{
		throw MatrixException(ERROR_IN_CREATING_MATRIX);
	}
}

/**
* A default destructor.
*/
template<class T>
Matrix<T>::~Matrix()
{}

/**
* Assignment operator
* @param matrix - the receives matrix.
* @return this after it was assigned by the receives matrix.
*/
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

/**
* + operator. Receives a matrix and returns a new matrix that represents the addition of
* the original matrix and the received matrix.
* @param matrix - the receives matrix.
* @return a new matrix that is the addition of the original and the receives matrix.
*/
template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &rhs) const
{
	if ((_rows != rhs._rows) || (_cols != rhs._cols))
	{
		throw MatrixException(ERROR_IN_PERFORMING_ADD_OPERATOR);
	}

	return _isParallel
		   ? _addParallel(*this, rhs)
		   : _addSerial(*this, rhs);
}

/**
* performs a regular add.
* @param lhs left matrix
* @param rhs right matrix
* @return the new matrix after the progress
*/
template<class T>
Matrix<T> Matrix<T>::_addSerial(const Matrix<T> &lhs, const Matrix<T> &rhs)
{
	Matrix<T> matrix(lhs._rows, lhs._cols);

	for (unsigned int i = 0; i < lhs._rows * lhs._cols; i++)
	{
		matrix._data[i] = lhs._data[i] + rhs._data[i];
	}
	return matrix;
}

/**
 * An object function. Receives three matrices in the constructor and has a () operator
 * which performs a single operation - row + row
 * @tparam T the type of the matrices elements
 */
template<typename T>
class VectorAdd
{
public:

	/**
	 * Receives three matrices - left, right and the new matrix to update
	 * @param lhs left matrix
	 * @param rhs right matrix
	 * @param result new matrix after the progress
	 */
	VectorAdd(const Matrix<T> &lhs, const Matrix<T> &rhs, Matrix<T> &result) :
			_lhs(lhs),
			_rhs(rhs),
			_result(result)
	{}

	/**
	 * performs a single operation for each thread - row + row
	 * @param i row
	 */
	void operator()(unsigned int i)
	{
		for (unsigned int j = 0; j < _lhs.cols(); j++)
		{
			_result(i, j) = _lhs(i, j) + _rhs(i, j);
		}
	}

private:

	const Matrix<T> &_lhs;
	const Matrix<T> &_rhs;
	Matrix<T> &_result;
};

/**
* Performs a parallel programming add
* @param lhs left matrix
* @param rhs right matrix
* @return the new matrix after the progress
*/
template<class T>
Matrix<T> Matrix<T>::_addParallel(const Matrix<T> &lhs, const Matrix<T> &rhs)
{
	Matrix<T> result(lhs._rows, rhs._cols);

	VectorAdd<T> adder(lhs, rhs, result);

	std::vector<std::thread> threads(lhs._rows);

	for (unsigned int i = 0; i < lhs._rows; i++)
	{
		threads[i] = std::thread(std::ref(adder), i);
	}

	for (auto &t : threads)
	{
		t.join();
	}
	return result;
}

/**
* - operator. Receives a matrix and returns a new matrix that represents the reduction of
* the original matrix and the received matrix.
* @param matrix - the receives matrix.
* @return a new matrix that is the reduction of the original and the receives matrix.
*/
template<class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &rhs) const
{
	if (_rows != rhs._rows || _cols != rhs._cols)
	{
		throw MatrixException(ERROR_IN_PERFORMING_REDUCTION_OPERATOR);
	}

	Matrix<T> matrix(_rows, _cols);

	for (unsigned int i = 0; i < _rows * _cols; i++)
	{
		matrix._data[i] = _data[i] - rhs._data[i];
	}

	return matrix;
}

/**
* operator. Receives a matrix and returns a new matrix that represents the multiplication of
* the original matrix and the receives matrix that the operator has received.
* @param matrix - the receives matrix.
* @return a new matrix that is the product of the original matrix and the receives matrix.
*/
template<class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &rhs) const
{
	if (_cols != rhs._rows)
	{
		throw MatrixException(ERROR_IN_PERFORMING_MULT_OPERATOR);
	}

	return _isParallel
		   ? _multParallel(*this, rhs)
		   : _multSerial(*this, rhs);
}

/**
* performs a regular multification.
* @param lhs left matrix
* @param rhs right matrix
*/
template<class T>
Matrix<T> Matrix<T>::_multSerial(const Matrix<T> &lhs, const Matrix<T> &rhs)
{
	Matrix<T> result(lhs._rows, rhs._cols);

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

/**
 * An object function. Receives three matrices in the constructor and has a () operator
 * which performs a single operation - row * col
 * @tparam T the type of the matrices elements
 */
template<typename T>
class VectorMultiplier
{
public:

	/**
	 * Receives three matrices - left, right and the new matrix to update
	 * @param lhs left matrix
	 * @param rhs right matrix
	 * @param result new matrix after the progress
	 */
	VectorMultiplier(const Matrix<T> &lhs, const Matrix<T> &rhs, Matrix<T> &result) :
			_lhs(lhs),
			_rhs(rhs),
			_result(result)
	{}

	/**
	 * performs a single operation for each thread - row * col
	 * @param i row
	 */
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
	Matrix<T> &_result;
};

/**
* Performs a parallel programming multification
* @param lhs left matrix
* @param rhs right matrix
* @return the new matrix after the progress
*/
template<class T>
Matrix<T> Matrix<T>::_multParallel(const Matrix<T> &lhs, const Matrix<T> &rhs)
{
	Matrix<T> result(lhs._rows, rhs._cols);

	VectorMultiplier<T> multiplier(lhs, rhs, result);

	std::vector<std::thread> threads(lhs._rows);

	for (unsigned int i = 0; i < lhs._rows; i++)
	{
		threads[i] = std::thread(std::ref(multiplier), i);
	}

	for (auto &t : threads)
	{
		t.join();
	}

	return result;
}

/**
* This operator receives a matrix and checks whether it is the same as the received matrix.
* @param matrix - the matrix we are comparing to.
* @return true if the matrix is equals to this, false otherwise.
*/
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

/**
* This operator receives a matrix and checks whether it is not the same as the received matrix.
* @param matrix - the matrix we are comparing to.
* @return true if the matrix is not the same as - this, false otherwise.
*/
template<class T>
bool Matrix<T>::operator!=(const Matrix<T> &rhs) const
{
	return !operator==(rhs);
}

/**
* @return true iff this is quadratic, false otherwise.
*/
template<class T>
bool Matrix<T>::isSquareMatrix() const
{
	return _rows == _cols;
}

/**
* trans() operator . Returns a new matrix that is the trans of - this.
* @return a new matrix that is the trans of - this.
*/
template<class T>
Matrix<T> Matrix<T>::trans() const
{
	if (!isSquareMatrix())
	{
		throw MatrixException(ERROR_IN_PERFORMING_TRANS);
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

/**
 * Specialization of trans() for type Complex.
 * @return the new matrix after performing the trans operation.
 */
template<>
inline Matrix<Complex> Matrix<Complex>::trans() const
{
	if (!isSquareMatrix())
	{
		throw MatrixException(ERROR_IN_PERFORMING_TRANS);
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

/**
* This friend operator receives a matrix and an output stream and prints the matrix.
* @tparam P - a class type.
* @param os - output stream.
* @param print - matrix to print.
* @return the received output stream.
*/
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

/**
* This operator receives a row and a column and returns the const element located in the [row,
* col] int this matrix.
* @param row - The row of the element.
* @param col - The column of the element.
* @return the element which is located in the given row and column
*/
template<class T>
T &Matrix<T>::operator()(unsigned int row, unsigned int col)
{
	if ((row > _rows) || (row < 0))
	{
		throw std::out_of_range(INDEX_OUT_OF_RANGE_OF_ROWS);
	}
	if ((col > _cols) || (col < 0))
	{
		throw std::out_of_range(INDEX_OUT_OF_RANGE_OF_COLS);
	}
	return _data[row * _cols + col];
}

/**
* This operator receives a row and a column and returns a const element located in the [row,col]
* int this matrix.
* @param row - The row of the element.
* @param col - The column of the element.
* @return the element which is located in the given row and column
*/
template<class T>
const T &Matrix<T>::operator()(unsigned int row, unsigned int col) const
{
	if ((row > _rows) || (row < 0))
	{
		throw std::out_of_range(INDEX_OUT_OF_RANGE_OF_ROWS);
	}
	if ((col > _cols) || (col < 0))
	{
		throw std::out_of_range(INDEX_OUT_OF_RANGE_OF_COLS);
	}
	return _data[row * _cols + col];
}

/**
* This method receives a bool value which determined the performance of the ", * operations.
* If is is true, the multiplication and plus operations of all the matrices will be performed
* by parallel programming.
* Otherwise (false value) will return the class to its default behavior (withous using
* parallel programming).
* @param parallel - a bool value which determines the mode of the operations +, *
* true if parallel, false if non-parallel.
*/
template<class T>
void Matrix<T>::setParallel(bool isParallel)
{
	if (isParallel != _isParallel)
	{
		std::string mode = isParallel ? PARALLEL_MODE : NON_PARALLEL_MODE;
		std::cout << CHANGE_MODE << mode << CHANGE_MODE_PART_2 << std::endl;
		_isParallel = isParallel;
	}
}


