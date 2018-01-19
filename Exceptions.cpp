#include "Exceptions.h"

Exceptions::Exceptions(const std::string &string) :
		_string(string)
{}

const char *Exceptions::what() const noexcept
{
	return _string.c_str();
}
