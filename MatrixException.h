#pragma once

#include <exception>

class MatrixException : public std::exception
{
public:
	virtual const char * what() const throw();
	virtual ~MatrixException();
};