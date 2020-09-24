/*
 * Solver_main.cpp
 *
 *  Created on: Jan 28, 2020
 *      Author: Nicolas
 */
/*
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




int main(){
//    int rows = 4;
//    int cols = 4;
//
////    double* b = new double[rows];
////	double* x = new double[rows];
//	std::vector<double> b(rows);
//	std::vector<double> x(rows);
//	b[0] = 8;
//	b[1] = 6;
//	b[2] = 4;
//	b[3] = 2;
//
//	for (int i = 0; i < rows; i++)
//		x[i] = 0;
//
//	//lets fill our matrix
//	auto *values = new double [rows * cols];
////	std::vector<double> values(rows * cols);
//		values[0] = 1.0;
//		values[1] = 2.0;
//		values[2] = 3.0;
//		values[3] = 4.0;
//		values[4] = 2.0;
//		values[5] = 3.0;
//		values[6] = 4.0;
//		values[7] = 6.0;
//		values[8] = 3.0;
//		values[9] = 4.0;
//		values[10] = 2.0;
//		values[11] = 5.0;
//		values[12] = 4.0;
//		values[13] = 6.0;
//		values[14] = 5.0;
//		values[15] = 7.0;
//
//	    double* LT_values = new double[rows*cols];
//
//	    for (int i = 0; i < rows * cols; i++)
//	    	LT_values[i] = 0;
//
//	    int A[] = {1,2,3,4,2,3,4,6,3,4,2,5,4,6,5,7};
//	    auto *dense_mat = new Matrix<double>(rows, cols, true); // Now we need to go and fill our matrix. Let's just fill it with A
//
//	    for (int i = 0; i < rows * cols; i++)
//	    {dense_mat->values[i] = A[i];}
//
//	    //Initialise L matrix
//	    auto *L_mat = new Matrix<double>(rows, cols, true); // Now we need to go and fill our matrix. Let's just fill it with A
//
//	    for (int i = 0; i < rows * cols; i++)
//	    {L_mat->values[i] = 0;}
//
//	auto* LT = new Matrix<double>(rows, cols, LT_values);
//	auto *mat = new Matrix<double>(rows, cols, A);
//
//	dense_mat->LU(*L_mat, b, x);
////    mat->LU(*dense_mat, b, x);
//    mat->printVector(x);
//    mat->Gauss_Elim(b, x);
//    mat->backsubstitution(b, x);
//
//    std::cout << "Reduced to upper triangular matrix: ";
//    std::cout << "Vector b: ";
//    mat->printVector(b);
//    std::cout << "Vector x: ";
//    mat->printVector(x);
//
//	b[0] = 8.0;
//	b[1] = 6.0;
//	b[2] = 4.0;
//	b[3] = 2.0;
//	for (int i = 0; i < rows; i++)
//		x[i] = 2.0;
//
////    std::cout << "Conjugate Gradient x: ";
////    mat->Conj_Grad(b, x);
////    mat->printVector(x);
////
////    std::cout << std::endl;
//
//
//    rows = 3;
//    cols = 3;
//    double* value = new double[rows*cols];
////    std::vector<double> value (rows * cols);
//	value[0] = 25.0;
//	value[1] = 15.0;
//	value[2] = 5.0;
//	value[3] = 15.0;
//	value[4] = 18.0;
//	value[5] = 0.0;
//	value[6] = 5.0;
//	value[7] = 0.0;
//	value[8] = 11.0;
//	auto* matrix = new Matrix<double>(rows, cols, value);
//
////	auto* d = new double[rows];
////	auto* xx = new double[rows];
//
//	std::vector<double> d(rows);
//	std::vector<double> xx(rows);
//	d[0] = 120.0;
//	d[1] = 84.0;
//	d[2] = 65.0;
//	for (int i = 0; i < rows; i++)
//		xx[i] = 2.0;
//
//	std::cout << std::endl << "Conjugate Gradient x: ";
//	matrix->Conj_Grad(d, xx);
//	matrix->printVector(xx);
//	std::cout << std::endl;
//
//
//	std::vector<double> y (rows);
//	for (int i = 0; i < rows; i++)
//		y[i] = 2.0;
//    auto* L_values = new double[rows * cols];
////    matrix->Cholesky(L_values);
//    auto* CL = new Matrix<double>(rows, cols, L_values);
//    CL->printMatrix();
//
////    double* LT_values = new double[rows*cols];
//    for (int i = 0; i < rows * cols; i++)
//    	LT_values[i] = 0;
//	for (int i = 0; i < rows; i++)
//		xx[i] = 2.0;
////    auto* LT = new Matrix<double>(rows, cols, LT_values);
//
////    CL->Transposer(*LT);
////    LT->printMatrix();
////    CL->Forward_Substitution(d, y);
////    LT->backsubstitution(y, xx);
////    std::cout << "xx : " << std::endl;
////    CL->printVector(xx);
////    matrix->Cholesky(*LT, d, xx);
////    matrix->printVector(xx);
//
//
//
//
//    ////////// SPARSE STUFF
//
//
////    int rows_sparse = 10;
////      int cols_sparse = 10;
////      double sparsity = 2/double(rows_sparse);
////      int nnzs = sparsity * rows_sparse * rows_sparse;
////
////      auto *test_csr_mat = new CSRMatrix<double>(rows_sparse, cols_sparse, nnzs, true);
////
////      test_csr_mat->fill_random_sparse();
////      std::cout << "Pseudorandomly-generated diagonally-dominant sparse matrix:" << std::endl << std::endl ;
////      test_csr_mat->printMatrix();
//
//    int rows_sparse = 4;
//    int cols_sparse = 4;
////    double sparsity = 2/double(rows_sparse);
//    int nnzs = 7;
//
//    auto *test_csr_mat = new CSRMatrix<double>(rows_sparse, cols_sparse, nnzs, true);
//
//    test_csr_mat->values[0] = 25;
//    test_csr_mat->values[1] = 1;
//    test_csr_mat->values[2] = 18;
//    test_csr_mat->values[3] = 1;
//    test_csr_mat->values[4] = 11;
//    test_csr_mat->values[5] = 1;
//    test_csr_mat->values[6] = 3;
//
//    test_csr_mat->col_index[0] = 0;
//    test_csr_mat->col_index[1] = 2;
//    test_csr_mat->col_index[2] = 1;
//    test_csr_mat->col_index[3] = 0;
//    test_csr_mat->col_index[4] = 2;
//    test_csr_mat->col_index[5] = 0;
//    test_csr_mat->col_index[6] = 3;
//
//    test_csr_mat->row_position[0] = 0;
//    test_csr_mat->row_position[1] = 2;
//    test_csr_mat->row_position[2] = 3;
//    test_csr_mat->row_position[3] = 5;
//    test_csr_mat->row_position[4] = 7;
//
//	b[0] = 54.0;
//	b[1] = 54.0;
//	b[2] = 46.0;
//	b[3] = 3.0;
//	for (int i = 0; i < rows_sparse; i++)
//		x[i] = 2.0;
//
//    test_csr_mat->Conj_Grad(b, x);
//    test_csr_mat->printVector(x);
//    test_csr_mat->printMatrix();
////	auto* x_coord = new int[nnzs];
//    std::vector<double> x_coord(nnzs);
//
//	for (int i = 0; i < rows + 1; i++){
//		for (int val_index = test_csr_mat->row_position[i]; val_index < test_csr_mat->row_position[i + 1]; val_index++){
//			x_coord[val_index] = i; // gives x coordinate of the value
//			std::cout << x_coord[val_index] << std::endl;
//		}
//	}
//
//
//
//////////////////////////////////////
//
////	nnzs = 4;
////    auto* L_sparse = new CSRMatrix<double>(rows_sparse, cols_sparse, nnzs, true);
////    L_sparse->row_position[0] = 0;
////    for (int i = 0; i < nnzs; i++){
////    	L_sparse->values[i] = 1;
////    	L_sparse->row_position[i] = i;
////    	L_sparse->col_index[i] = i;
////    }
////    L_sparse->row_position[nnzs] = nnzs;
////    L_sparse->printMatrix();
//
//    std::cout << "Sparse Cholesky: " << std::endl;
////    test_csr_mat->Cholesky();
////    L_sparse->printMatrix();
//
//
//
//
//    ////////SPD DENSE TEST


//    b.clear();
//    x.clear();
//    auto* x1 = new double[cols];
//    auto* b1 = new double[cols];
//    auto* blank = new double[cols];




//    SPDdense->fill_spd_dense();
//    SPDdense->printMatrix();
//    SPDdense->Conj_Grad(b1, x1);
//    std::cout << "x1 :" << std::endl;
//    SPDdense->printVector(x1);
//    SPDdense->matVecMult(x1, blank);
//    SPDdense->printVector(blank);

/////////// SPD CG SPARSE TeST

//    for (int i = 0; i < cols; i++)
//    {
////      b1[i] = rand() % 10;
//      x[i] = 0;
////      blank[i] = 0;
//    }
//    CSRMatrix<double> *SPDsparse = generate_random_spd_sparse(rows, cols);
//    SPDsparse->printMatrix();
//    SPDsparse->Conj_Grad(b1, x1);
//    SPDsparse->matVecMult(x1, blank);
//    SPDsparse->printVector(blank);



//    std::cout << "///////////////////////////////////////////////////////////////" << std::endl;
//    test_csr_mat->Cholesky();

      std::cout << "Matrix Library !!!";
      std::cout << "There are a variety of solvers here and this library can demonstrate their use,";
      std::cout << "you can input a square matrix row size and see how they perform.";





      int rows = 15;
      int cols = 15;

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


*/
