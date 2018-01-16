#include "MatrixException.h"

const char* MatrixException::what() const
{
	return "invalid matrix cordinates";
}

MatrixException::~MatrixException()
{
	delete MatrixException;
}
