#include "MatrixException.h"

MatrixException::MatrixException(const std::string &string) :
		_string(string)
{}

const char *MatrixException::what() const noexcept
{
	return _string.c_str();
}
