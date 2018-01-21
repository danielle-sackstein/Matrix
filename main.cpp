

#include "Matrix.hpp"

using std::cout;
using std::endl;

int main()
{
	//Matrix<int>::setParallel(true);
	cout << "\n*Tests number 3, 5, 6 compares between doubles, according to Yoni the Metargel, problems may occur, if you failed this test look your printed result because you might be alright, this process result doesn't always have consistent behaviour*\n" << endl;

	// test 1
	cout << "Test 1: checks complex addition: ";
	Matrix<Complex> comp1(5, 3, {9, 8, 7, 1, 2, 3, 4, 4, 4, 3, 3, 3, 5, 4, 3});
	Matrix<Complex> comp2(5, 3, {6, 7, 7, 4, 4, 2, 8, 8, 1, 9, 9, 0, 5, 5, 6});
	Matrix<Complex> expectedResult1(5, 3, {15, 15, 14, 5, 6, 5, 12, 12, 5, 12, 12, 3, 10, 9, 9});
	Matrix<Complex> realResult1 = comp1 + comp2;

	if(expectedResult1 != realResult1)
	{
		cout << "Failed" << endl;
		cout << "Expected: " << expectedResult1 << endl << "Result: " << realResult1 << endl;

	}else{
		cout << "Passed" << endl;
	}

	// test 2
	cout << "\nTest 2: checks complex trans: ";
	Matrix<Complex> cmplx(4, 4, {1, 8, 8, 8, 4, 2, 4, 4, 3, 3, 5, 3, 9, 9, 9, 0});
	Matrix<Complex> expectedResult2(4, 4, {1, 4, 3, 9, 8, 2, 3, 9, 8, 4, 5, 9, 8, 4, 3, 0});
	Matrix<Complex> realResult2 = cmplx.trans();
	if(expectedResult2 != realResult2)
	{
		cout << "Failed" << endl;
		cout << "Expected: " << expectedResult2 << endl << "Result: " << realResult2 << endl;

	}else{
		cout << "Passed" << endl;
	}

	// test 3
	cout << "\nTest 3: checks double multiplication: ";
	Matrix<double> doubleM1(3, 3, {9, 9, 9, 8, 8, 8, 7, 7, 7});
	Matrix<double> doubleM2(3, 3, {1, 0, 1, 2, 1, 2, 3, 2, 3});
	Matrix<double> expectedResult3(3, 3, {54, 27, 54, 48, 24, 48, 42, 21, 42});
	Matrix<double> realResult3 = doubleM1 * doubleM2;
	if(expectedResult3 != realResult3)
	{
		cout << "Failed" << endl;
		cout << "Expected: " << expectedResult3 << endl << "Result: " << realResult3 << endl;

	}else{
		cout << "Passed" << endl;
	}

	// test 4
	cout << "\nTest 4: checks complex multiplication: ";
	Matrix<Complex> c1(4, 4, {1, 0, 1, 0, 2, 1, 2, 1, 3, 2, 3, 2, 4, 3, 4, 3});
	Matrix<Complex> c2(4, 4, {9, 8, 9, 8, 8, 7, 8, 7, 6, 5, 6, 5, 5, 4, 5, 4});
	Matrix<Complex> expectedResult4(4, 4, {15, 13, 15, 13, 43, 37, 43, 37, 71, 61, 71, 61, 99, 85, 99, 85});
	Matrix<Complex> realResult4 = c1 * c2;
	if(expectedResult4 != realResult4)
	{
		cout << "Failed" << endl;
		cout << "Expected: " << expectedResult4 << endl << "Result: " << realResult4 << endl;
	}else{
		cout << "Passed" << endl;
	}

	// test 5
	cout << "\nTest 5: checks double multiplication: ";
	Matrix<double> d1(2, 2, {8, 8, 2, 2});
	Matrix<double> d2(2, 2, {6, 6, 3, 3});
	Matrix<double> expected5(2, 2, {72, 72, 18, 18});
	Matrix<double> realResult5 = d1*d2;
	if(expected5 != realResult5)
	{
		cout << "Failed" << endl;
		cout << "Expected: " << expected5 << endl << "Result: " << realResult5 << endl;
	}else{
		cout << "Passed" << endl;
	}

	// test 6
	cout << "\nTest 6: checks double addition: ";
	Matrix<double> db1(4, 2, {7, 7, 0, 0, 6, 6, 1.5, 2.5});
	Matrix<double> db2(4, 2, {8.8, 9.9, 7.7, 6.6, 5.5, 4.4, 3.3, 2.1});
	Matrix<double> expected6(4, 2, {15.8, 16.9, 7.7, 6.6, 11.5, 10.4, 4.8, 4.6});
	Matrix<double> realResult6 = db1 + db2;
	if(expected6 != realResult6)
	{
		cout << "Failed" << endl;
		cout << "Expected: \n" << expected6 << endl << "Result: \n" << realResult6 << endl;
	}else{
		cout << "Passed" << endl;
	}

	// test 7
	cout << "\nTest 7: checks if addition of different size matrices<int> fails correctly, your program suppose to throw std::exception inheritor in that case. result:\n";
	Matrix<int> im1(2, 3);
	Matrix<int> im2(4, 7);
	try
	{
		std::cerr << im1 + im2 << "\n you didn't throw something that inherits from std::exception!\n";
	}catch (std::exception e) {
		cout << "Exception thrown !! OK :)\n";
	}

	// test 8
	cout << "\nTest 8: checks access to out of matrix<int> bounds - this test checks if you throw std::exception inheritor. result:\n";
	Matrix<int> im3(4, 5, {2, 2, 3});
	try
	{
		std::cerr << im3(6, 6) << "\n you suppose to throw std::exception inheritor\n";
	}catch (std::exception e) {
		cout << "Exception thrown !! OK :)\n";
	}

	// test 9
	cout << "\nLast test 9 checks matrices<int> subtraction, which school sol doesn't allow, this is my result, let me know if you convinced im wrong (natalia):\n";
	Matrix<int> test1(3, 3, {1, 2, 3, 4, 5, 6, 7, 8, 9});
	Matrix<int> test2(3, 3, {9, 8, 7, 6, 5, 4, 3, 2, 1});
	Matrix<int> res(3, 3, {-8, -6, -4, -2, 0, 2, 4, 6, 8});
	cout << "My result is: \n" << res << endl;
	cout << "Your result is: \n" << test1 - test2 << endl;

	cout << "\nGOOD LUCK EVERYONE! :]\n";

	return 0;
}
