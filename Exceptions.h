#pragma once

#include <exception>
#include <string>

class Exceptions : public std::exception
{
public:
	Exceptions(const std::string &);

	const char *what() const noexcept override;

private :
	std::string _string;
};