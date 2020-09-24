/*
 * Main.cpp
 *
 *  Created on: Feb 2, 2020
 *      Author: Nicolas
 */

#include <iostream>
#include <vector>
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator
#include <sstream>
#include "Matrix.cpp"
#include "Matrix.h"
#include "CSRMatrix.cpp"
#include "CSRMatrix.h"
#include "Matrix_Tester.cpp"
#include "Matrix_Tester.h"

template <class T>
void solver(Matrix<T> M, std::vector<T> b, std::vector<T> x) {
	M.jacobi(b, x);
	M.printVector(b);
}

int main() {

    std::cout << "Matrix Library !!!";
    std::cout << "There are a variety of solvers here and this library can demonstrate their use," << std::endl;
    std::cout << "you can input a square matrix row size and see how they perform." << std::endl << std::endl;

    std::cout << "Here, we solve an arbitrary 20 by 20 dense diagonally dominant matrix with Jacobi: " << std::endl;



    auto* Mat = new Matrix<double>(20, 20, true);
    Mat->fill_random_dense();
	std::vector<double> x_arb(20);
	std::vector<double> b_arb(20);
	std::vector<double> comp(20);
	for (int i = 0; i < 20; i++){
		b_arb[i] = rand() % 10;
		x_arb[i] = 1;
	}
	std::cout << "Here is an arbitrary 20 by 20 dense matrix: " << std::endl;
	Mat->printMatrix();
	std::cout << "Here is an arbitrary b vector: " << std::endl;
	Mat->printVector(b_arb);
	std::cout << "Here is a guess x vector: " << std::endl;
	Mat->printVector(x_arb);
	std::cout << "Now, let's see the result" << std::endl;
    solver(*Mat, b_arb, x_arb);


    std::cout << "Are we sure that it is correct? We can find out by multiplying the";
    std::cout << " result with the matrix to see if we get the same b vector.";
	std::cout << std::endl << "Checking result...";
	Mat->matVecMult(x_arb, comp);

	std::cout << "We obtain this vector when checking the result: ";
    Mat->printVector(comp);
	for (int i = 0; i < 20; i++)
		if ( b_arb[i] - comp[i] > 1e-3 ){
			std::cout << "Wrong";
		}
	std::cout << "The solution was correct!" << std::endl;
	delete Mat;


	int rows = 15;
	int cols = 15;
	auto* Example_Mat = new Matrix<double>(rows, cols, true);
	std::cout << "Diagonally dominant random dense matrix." << endl;
	Example_Mat->fill_random_dense();
	Example_Mat->printMatrix();


	std::cout << "Diagonally dominant CSR matrix." << std::endl;
	double sparsity = 0.05;
	int nnzs = sparsity * rows * cols;
	CSRMatrix<double> *test_csr_mat = new CSRMatrix<double>(rows, cols, nnzs, true);
    auto* SPDsparse = test_csr_mat->generate_random_spd_sparse(rows, cols);

	test_csr_mat->printMatrix();


	int ml = 1;
	int mu = 2;
	std::cout << "This is a pseudo-randomly generated banded matrix with ml = " << ml << ", mu = " << mu << ", and bandwidth = " << mu + ml + 1 << ":" << std::endl;
	BandedMatrix<double> *test_banded_mat = new BandedMatrix<double>(rows, cols, ml, mu, true);

	test_banded_mat->fill_random_banded();
	test_banded_mat->printMatrix();





	///////// Test for all solvers with random matrices for each.
	rows = 15;
	cols = 15;

	std::vector<double> x(cols);
	std::vector<double> b(cols);

	for (int i = 0; i < cols; i++){
		b[i] = rand() % 10;
		x[i] = 1;
	}
	std::cout << std::endl << "b vector (for solving linear system Ax = b):" << std::endl;
	for (int i = 0; i < cols; i++)
		std::cout << b[i] << " ";

	auto* TEST = new Matrix_Tester<double>();


	///SPD Tester
	Matrix<double>* SPDdense = new Matrix<double>(rows, cols, true);
	SPDdense->fill_spd_dense();
	bool dense_is_spd = TEST->SPD_dense_test(*SPDdense);



	///Solver Testing
	TEST->Gauss_Elim_test(rows, cols, b);
	std::cout << std::endl<< std::endl << std::endl;
	TEST->CG_test(rows, cols, false, b);
	std::cout << std::endl<< std::endl << std::endl;
	TEST->CG_test(rows, cols, true, b);
	std::cout << std::endl<< std::endl << std::endl;
	TEST->jacobi_test(rows, cols, true, true, b);
	std::cout << std::endl<< std::endl << std::endl;
	TEST->Gauss_Seidel_test(rows, cols, true, true, b);
	std::cout << std::endl<< std::endl << std::endl;
	TEST->SOR_test(rows, cols, true, true, b);
	std::cout << std::endl<< std::endl << std::endl;
	TEST->Cholesky_test(rows, cols, b);



	return 0;
}
