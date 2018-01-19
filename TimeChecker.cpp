//#include <chrono>
//#include <stack>
//#include "Matrix.hpp"
//#include <Eigen/Dense>
//
//using Eigen::MatrixXd;
//
//std::stack<std::chrono::time_point<std::chrono::system_clock>> tictoc_stack;
//
///**
// * starts the timer
// */
//void tic()
//{
//	tictoc_stack.push(std::chrono::system_clock::now());
//}
//
///**
// * ends the timer
// * @param message the message to print.
// */
//void toc(std::string message)
//{
//	std::chrono::duration<double> elapsed_seconds =
//			std::chrono::system_clock::now() - tictoc_stack.top();
//	std::cout << message << elapsed_seconds.count() << std::endl;
//	tictoc_stack.pop();
//}
//
//
//void printLibraryTime(unsigned int n)
//{
//	std::vector<int> v(n * n, 1);
//
//	Matrix<int> m1(n, n, v);
//	Matrix<int> m2(n, n, v);
//
//	// operation +
//	tic();
//
//	Matrix<int> add = m1 + m2;
//
//	toc("matlib add ");
//
//	// operation *
//	tic();
//
//	Matrix<int> mul = m1 * m2;
//
//	toc("matlib mult ");
//}
//
///**
// * calculates and prints the time of the operations + and * on the matrices on the library
// * Matrix.hpp
// * @param n the size of the matrix
// */
//void printParallelTime(unsigned int n)
//{
//	std::vector<int> v(n * n, 1);
//
//	Matrix<int> m1(n, n, v);
//	Matrix<int> m2(n, n, v);
//
//	// operation +
//	tic();
//
//	Matrix<int> add = m1 + m2;
//
//	toc("matlibP add ");
//
//	// operation *
//	tic();
//
//	Matrix<int> mul = m1 * m2;
//
//	toc("matlibP mult ");
//}
//
///**
// * calculates and prints the time of the operations + and * on the matrices on the library Eigen
// * @param n the size of the matrix
// */
//void printEigenTime(unsigned int n)
//{
//	MatrixXd m1 = MatrixXd::Random(n, n);
//	MatrixXd m2 = MatrixXd::Random(n, n);
//
//	// operation +
//	tic();
//
//	MatrixXd add = m1 + m2;
//
//	toc("eigen add ");
//
//	// operation *
//	tic();
//
//	MatrixXd mul = m1 * m2;
//
//	toc("eigen mult ");
//}
//
///**
// * main method. Receives a size of a matrix and creates matrices, performing + and * operations.
// * @param argc number of args
// * @param argv argument which represents the size of the matrix.
// * @return
// */
//int main(int argc, char *argv[])
//{
//	unsigned int n;
//	sscanf(argv[1], "%d", &n);
//
//	std::cout << "size " << n << std::endl;
//
//	printEigenTime(n);
//	printLibraryTime(n);
//	Matrix<int>::setParallel(true);
//	printParallelTime(n);
//
//}
//
