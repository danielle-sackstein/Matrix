#include "MatrixException.h"

/**
 * MatrixException class. receives  a string which represents the message to print
 * @param string  represents the message to print
 */
MatrixException::MatrixException(const std::string &message) :
		_message(message)
{}

/**
 * @return returns the message to print
 */
const char *MatrixException::what() const noexcept
{
	return _message.c_str();
}