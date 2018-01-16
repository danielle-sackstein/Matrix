#include "MatrixException.h"

char const *MatrixException::what() const
{
	return "invalid matrix cordinates";
}

MatrixException::~MatrixException()
{
	delete MatrixException;
}
