#include "MatrixException.h"

MatrixException::MatrixException(const std::string &string) :
		_message(string)
{}

const char *MatrixException::what() const noexcept
{
	return _message.c_str();
}