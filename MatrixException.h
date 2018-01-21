#pragma once

#include <exception>
#include <string>

class MatrixException : public std::exception
{
public:
	MatrixException(const std::string &);

	const char *what() const noexcept override;

private :
	std::string _message;
};