/*
 * Matrix_Tester.cpp
 *
 *  Created on: Feb 2, 2020
 *      Author: Nicolas
 */

#include "Matrix_Tester.h"
#include <iostream>
#include <vector>
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator
#include <sstream>
#include "BandedMatrix.cpp"

template <class T>
void Matrix_Tester<T>::CG_test(int rows, int cols, bool sparsed, std::vector<T>& b_val){
	std::vector<T> x_val(cols, 0); // initialize guess
	std::vector<T> comp(cols, 0); // compares the result with the initial b_val

	if (!sparsed){ // use SPD dense generators here
	    Matrix<T>* SPDdense = new Matrix<T>(rows, cols, true);

		SPDdense->fill_spd_dense(); // make it SPD
		SPDdense->Conj_Grad(b_val, x_val);
		std::cout << std::endl << "The dense Conjugate Gradient solution x is: ";
		SPDdense->printVector(x_val);
		std::cout << std::endl << "Checking result...";
		SPDdense->matVecMult(x_val, comp);
		SPDdense->printVector(comp);
	    std::cout << std::endl << "With b: ";
	    SPDdense->printVector(b_val);

		// check if b and comp are close to the same
		for (int i = 0; i < cols; i++)
			if ( b_val[i] - comp[i] > 1e-3){
				std::cout << "Wrong";
			}
		std::cerr << "Correct"<< std::endl;
		delete SPDdense;
	}
	else{

		auto* sparse = new CSRMatrix<T>(rows, cols, 1, true);
	    auto* SPDsparse = sparse->generate_random_spd_sparse(rows, cols);
	    SPDsparse->Conj_Grad(b_val,x_val);
	    SPDsparse->matVecMult(x_val, comp);
		std::cout << std::endl << "The sparse Conjugate Gradien solution x is: ";
		SPDsparse->printVector(x_val);
		std::cout << std::endl << "Checking result...";
	    SPDsparse->printVector(comp);
	    std::cout << std::endl << "With b: ";
	    SPDsparse->printVector(b_val);

		for (int i = 0; i < cols; i++)
			if ( b_val[i] - comp[i] > 1e-3 ){
				std::cerr << "Wrong";
			}
		std::cout << "Correct" << std::endl;
		delete sparse;
		delete SPDsparse;
	}
}
template <class T>
void Matrix_Tester<T>::jacobi_test(int rows, int cols, bool sparsed, bool banded, std::vector<T>& b_val){

	std::vector<T> x_val(cols, 0); // initialize guess
	std::vector<T> comp(cols, 0); // compares the result with the initial b_val

	if(!sparsed && !banded) {
	    Matrix<T>* DDdense = new Matrix<T>(rows, cols, true);

		DDdense->fill_random_dense_diagdom(); // make it diagonally dominant
		DDdense->jacobi(x_val, b_val);
		DDdense->matVecMult(x_val, comp);
		std::cout << std::endl << "The dense Jacobi solution x is: ";
		DDdense->printVector(x_val);
		std::cout << std::endl << "Checking result...";
	    DDdense->printVector(comp);
	    std::cout << std::endl << "With b: ";
	    DDdense->printVector(b_val);
	    DDdense->printVector(comp);

		for (int i = 0; i < cols; i++)
			if ( b_val[i] - comp[i] > 1e-3 ){
				std::cout << "Wrong";
			}
		std::cout << "Correct" << std::endl;
		delete DDdense;

	}
	else{
		if (!banded){
			double sparsity;
			std::cout << "Choose amount of sparsity: ";
			std::cin >> sparsity;


			  // determines number of nonzeros to add to the sparse matrix. Note that since this method generates a matrix with a uniform number of nonzero elements in each row, it adjusts the user specified sparsity accordingly.

			int nnzs;
			while (sparsity <= 0 || sparsity >= 0.25){
				if (sparsity <= 0 || sparsity >= 1)
				    std::cerr << "Invalid sparsity, which must be between 0 and 1" << std::endl;
				if (1 > sparsity > 0.25)
					std::cerr << "This matrix is quite dense. Please enter a lower sparsity." << std::endl;
			  double nnzs_input = (rows * cols) * sparsity;
			  nnzs = round(nnzs_input / static_cast<double>(rows)) * static_cast<double>(rows);

			  if (nnzs < rows)
			  {
				  std::cout << "Not enough elements to fill the diagonal entries. Please enter a greater sparsity." << std::endl;
			  }
			  std::cin >> sparsity;

			  }
			auto* sparse = new CSRMatrix<T>(rows, cols, nnzs, true);
			auto* diagdom_sparse = sparse->generate_random_diagdom_sparse(rows, cols, sparsity);
			diagdom_sparse->jacobi(x_val, b_val);
			diagdom_sparse->matVecMult(x_val, comp);
			std::cout << std::endl << "The sparse Jacobi solution x is: ";
			diagdom_sparse->printVector(x_val);
			std::cout << std::endl << "Checking result...";
		    diagdom_sparse->printVector(comp);
		    std::cout << std::endl << "With b: ";
		    diagdom_sparse->printVector(b_val);
		    diagdom_sparse->printVector(comp);
			for (int i = 0; i < cols; i++)
				if ( b_val[i] - comp[i] > 1e-3 ){
					std::cout << "Wrong";
				}
			std::cout << "Correct" << std::endl;
			delete sparse;
			delete diagdom_sparse;
		}
		else {
			int ml = 1;
			int mu = 2;
			BandedMatrix<double> *test_banded_mat = new BandedMatrix<double>(rows, cols, ml, mu, true);
			test_banded_mat->fill_random_banded();
			cout << "This is a pseudo-randomly generated banded matrix with ml = " << ml << ", mu = " << mu << ", and bandwidth = " << mu + ml + 1 << ":" << endl;
			test_banded_mat->jacobi(x_val, b_val);
			test_banded_mat->matVecMult(x_val, comp);
			std::cout << std::endl << "The banded Jacobi solution x is: ";
			test_banded_mat->printVector(x_val);
			std::cout << std::endl << "Checking result...";
			test_banded_mat->printVector(comp);
			std::cout << std::endl << "With b: ";
			test_banded_mat->printVector(b_val);
			test_banded_mat->printVector(comp);
			for (int i = 0; i < cols; i++)
				if ( b_val[i] - comp[i] > 1e-3 ){
					std::cout << "Wrong";
				}
			std::cout << "Correct" << std::endl;
			delete test_banded_mat;
//			fill(x.begin(), x.end(), 1);
//			fill(output.begin(), output.end(), 0);
		}
	}
}


template <class T>
void Matrix_Tester<T>::Gauss_Seidel_test(int rows, int cols, bool sparsed, bool banded, std::vector<T>& b_val){
	std::vector<T> x_val(cols, 0); // initialize guess
	std::vector<T> comp(cols, 0); // compares the result with the initial b_val

	if(!sparsed && !banded) {
		Matrix<T>* DDdense = new Matrix<T>(rows, cols, true);

		DDdense->fill_random_dense_diagdom(); // make it diagonally dominant
		DDdense->Gauss_Seidel(x_val, b_val);
		DDdense->matVecMult(x_val, comp);
		std::cout << std::endl << "The dense GS solution x is: ";
		DDdense->printVector(x_val);
		std::cout << std::endl << "Checking result...";
		DDdense->printVector(comp);
		std::cout << std::endl << "With b: ";
		DDdense->printVector(b_val);
		DDdense->printVector(comp);

		for (int i = 0; i < cols; i++)
			if ( b_val[i] - comp[i] > 1e-3 ){
				std::cout << "Wrong";
			}
		std::cout << "Correct" << std::endl;
		delete DDdense;

	}
	else{
		if (!banded){ // SPARSE
			int nnzs;
			double sparsity = 0.1;
			double nnzs_input = (rows * cols) * sparsity;
			nnzs = round(nnzs_input / static_cast<double>(rows)) * static_cast<double>(rows);

			auto* sparse = new CSRMatrix<T>(rows, cols, nnzs, true);
			auto* diagdom_sparse = sparse->generate_random_diagdom_sparse(rows, cols, sparsity);
			diagdom_sparse->Gauss_Seidel(x_val, b_val);
			diagdom_sparse->matVecMult(x_val, comp);
			std::cout << std::endl << "The sparse GS solution x is: ";
			diagdom_sparse->printVector(x_val);
			std::cout << std::endl << "Checking result...";
			diagdom_sparse->printVector(comp);
			std::cout << std::endl << "With b: ";
			diagdom_sparse->printVector(b_val);
			diagdom_sparse->printVector(comp);
			for (int i = 0; i < cols; i++)
				if ( b_val[i] - comp[i] > 1e-3 ){
					std::cout << "Wrong";
				}
			std::cout << "Correct" << std::endl;
			delete sparse;
			delete diagdom_sparse;
		}
		else{ //BANDED
			int ml = 1;
			int mu = 2;
			BandedMatrix<double> *DDBanded = new BandedMatrix<double>(rows, cols, ml, mu, true);

			DDBanded->fill_random_banded(); // make it diagonally dominant
			DDBanded->Gauss_Seidel(x_val, b_val);
			DDBanded->matVecMult(x_val, comp);
			std::cout << std::endl << "The banded GS solution x is: ";
			DDBanded->printVector(x_val);
			std::cout << std::endl << "Checking result...";
			DDBanded->printVector(comp);
			std::cout << std::endl << "With b: ";
			DDBanded->printVector(b_val);
			DDBanded->printVector(comp);

			for (int i = 0; i < cols; i++)
				if ( b_val[i] - comp[i] > 1e-3 ){
					std::cout << "Wrong";
				}
			std::cout << "Correct" << std::endl;
			delete DDBanded;
		}
	}
}



template <class T>
void Matrix_Tester<T>::SOR_test(int rows, int cols, bool sparsed, bool banded, std::vector<T>& b_val){
	std::vector<T> x_val(cols, 0); // initialize guess
	std::vector<T> comp(cols, 0); // compares the result with the initial b_val

	double w = 0.9;

	if(!sparsed && !banded) {
		Matrix<T>* DDdense = new Matrix<T>(rows, cols, true);

		DDdense->fill_random_dense_diagdom(); // make it diagonally dominant
		DDdense->sor(x_val, b_val, w);
		DDdense->matVecMult(x_val, comp);
		std::cout << std::endl << "The dense SOR solution x is: ";
		DDdense->printVector(x_val);
		std::cout << std::endl << "Checking result...";
		DDdense->printVector(comp);
		std::cout << std::endl << "With b: ";
		DDdense->printVector(b_val);
		DDdense->printVector(comp);

		for (int i = 0; i < cols; i++)
			if ( b_val[i] - comp[i] > 1e-3 ){
				std::cout << "Wrong";
			}
		std::cout << "Correct" << std::endl;
		delete DDdense;

	}
	else{
		if (!banded){ // SPARSE
			int nnzs;
			double sparsity = 0.1;
			double nnzs_input = (rows * cols) * sparsity;
			nnzs = round(nnzs_input / static_cast<double>(rows)) * static_cast<double>(rows);

			auto* sparse = new CSRMatrix<T>(rows, cols, nnzs, true);
			auto* diagdom_sparse = sparse->generate_random_diagdom_sparse(rows, cols, sparsity);
			diagdom_sparse->sor(x_val, b_val, w);
			diagdom_sparse->matVecMult(x_val, comp);
			std::cout << std::endl << "The sparse SOR solution x is: ";
			diagdom_sparse->printVector(x_val);
			std::cout << std::endl << "Checking result...";
			diagdom_sparse->printVector(comp);
			std::cout << std::endl << "With b: ";
			diagdom_sparse->printVector(b_val);
			diagdom_sparse->printVector(comp);
			for (int i = 0; i < cols; i++)
				if ( b_val[i] - comp[i] > 1e-3 ){
					std::cout << "Wrong";
				}
			std::cout << "Correct" << std::endl;
			delete sparse;
			delete diagdom_sparse;
		}
		else{ //BANDED
			int ml = 1;
			int mu = 2;
			BandedMatrix<double> *DDBanded = new BandedMatrix<double>(rows, cols, ml, mu, true);

			DDBanded->fill_random_banded(); // make it diagonally dominant
			DDBanded->sor(x_val, b_val, w);
			DDBanded->matVecMult(x_val, comp);
			std::cout << std::endl << "The banded SOR solution x is: ";
			DDBanded->printVector(x_val);
			std::cout << std::endl << "Checking result...";
			DDBanded->printVector(comp);
			std::cout << std::endl << "With b: ";
			DDBanded->printVector(b_val);
			DDBanded->printVector(comp);

			for (int i = 0; i < cols; i++)
				if ( b_val[i] - comp[i] > 1e-3 ){
					std::cout << "Wrong";
				}
			std::cout << "Correct" << std::endl;
			delete DDBanded;
		}
	}
}




template <class T>
bool Matrix_Tester<T>::SPD_dense_test(Matrix<T>& candidate_spd_dense_mat){
// If a matrix A is SPD, then for any nonzero vector z and its
// transpose z**T, z**T A z = positive scalar

	for (int i = 0; i < 100; i++)
	{
		std::vector<double> z(candidate_spd_dense_mat.cols, 0);
		std::vector<double> y(candidate_spd_dense_mat.cols, 0);
		double scalar = 0.0;

		// Generates a random nonzero vector z
		for (int i = 0; i < candidate_spd_dense_mat.cols; i++)
		{
			z[i] = rand() % 20 - 10;
		}

		// z**T A yields another row vector, y
		candidate_spd_dense_mat.matRowVecMult(z, y);

		for (int i = 0; i < z.size(); i++)
		{
			scalar += y[i] * z[i];
		}

		if (scalar < 0)
		{
			std::cout << "Not symmetric positive definite" << std::endl;
			return false;
		}
	}

	std::cout << "Dense Symmetric positive definite" << std::endl;
	return true;
}

template <class T>
void Matrix_Tester<T>::Cholesky_test(int rows, int cols, std::vector<T>& b_val){
	std::vector<T> x_val(cols, 0); // initialize guess
	std::vector<T> comp(cols, 0); // compares the result with the initial b_val
    T* LT_values = new T[rows*cols];
    for (int i = 0; i < rows * cols; i++)
    	LT_values[i] = 0;
	auto* L = new Matrix<T>(rows, cols, LT_values);

	Matrix<T>* SPDdense = new Matrix<T>(rows, cols, true);

	SPDdense->fill_spd_dense(); // make it diagonally dominant
	SPDdense->Cholesky(*L, b_val, x_val);
	SPDdense->matVecMult(x_val, comp);
	std::cout << std::endl << "The dense Cholesky solution x is: ";
	SPDdense->printVector(x_val);
	std::cout << std::endl << "Checking result...";
	SPDdense->printVector(comp);
	std::cout << std::endl << "With b: ";
	SPDdense->printVector(b_val);
	SPDdense->printVector(comp);

	for (int i = 0; i < cols; i++)
		if ( b_val[i] - comp[i] > 1e-3 ){
			std::cout << "Wrong";
		}
	std::cout << "Correct" << std::endl;
	delete SPDdense;
}

template <class T>
void Matrix_Tester<T>::Gauss_Elim_test(int rows, int cols, std::vector<T>& b_val){
	std::vector<T> x_val(cols, 0); // initialize guess
	std::vector<T> comp(cols, 0); // compares the result with the initial b_val
	Matrix<T>* dense = new Matrix<T>(rows, cols, true);

	dense->fill_spd_dense(); // make it diagonally dominant
	dense->Gauss_Elim(b_val, x_val);
	dense->backsubstitution(b_val, x_val);
	dense->matVecMult(x_val, comp);
	std::cout << std::endl << "The dense Gaussian elimination solution x is: ";
	dense->printVector(x_val);
	std::cout << std::endl << "Checking result...";
	dense->printVector(comp);
	std::cout << std::endl << "With b: ";
	dense->printVector(b_val);
	dense->printVector(comp);

	for (int i = 0; i < cols; i++)
		if ( b_val[i] - comp[i] > 1e-3 ){
			std::cout << "Wrong";
		}
	std::cout << "Correct" << std::endl;
	delete dense;

}

//
//template <class T>
//void Matrix_Tester<T>::LU_test(int rows, int cols, std::vector<T>& b_val){
//
//}
//
















