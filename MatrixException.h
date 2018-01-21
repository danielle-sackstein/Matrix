#pragma once

#include <exception>
#include <string>


class MatrixException : public std::exception
{
public:

	/**
	 * MatrixException class. receives  a string which represents the message to print
	 * @param string  represents the message to print
	 */
	MatrixException(const std::string &);

	/**
 	* @return returns the message to print
 	*/
	const char *what() const noexcept override;

private :
	std::string _message;
};