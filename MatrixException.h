#pragma once

#include <exception>

class MatrixException : public std::exception
{
public:
	virtual char const* what() const throw();
	virtual ~MatrixException();
};