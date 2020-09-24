/*
 * CSRMatrix.h
 *
 *  Created on: Jan 30, 2020
 *      Author: Nicolas
 */

#ifndef CSRMATRIX_H_
#define CSRMATRIX_H_

#pragma once
#include "Matrix.h"


template <class T>
class CSRMatrix: public Matrix<T>
{
public:

   // constructor where we want to preallocate ourselves
   CSRMatrix(int rows, int cols, int nnzs, bool preallocate);
   // constructor where we already have allocated memory outside
   CSRMatrix(int rows, int cols, int nnzs, T *values_ptr, int *row_position, int *col_index);
   // destructor
   ~CSRMatrix();

   // Print out the values in our matrix
   virtual void printMatrix();
   virtual void printVector(std::vector<double>& v);

   CSRMatrix<T>* transpose();

   static bool Sorter(const std::pair<int, int>& left, const std::pair<int, int>& right);

   void fill_random_sparse();

   // Perform some operations with our matrix
   CSRMatrix<T>* matMatMult(CSRMatrix<T>& B);
   CSRMatrix<T>* generate_random_spd_sparse(int rows, int cols);
   CSRMatrix<T>* generate_random_diagdom_sparse(int rows, int cols, double sparsity);
   virtual void matVecMult(std::vector<double>& input, std::vector<double>& output);
   virtual void matRowVecMult(std::vector<double>& input, std::vector<double>& output);

   // Explicitly using the C++11 nullptr here
   int *row_position = nullptr;
   int *col_index = nullptr;

   // How many non-zero entries we have in the matrix
   int nnzs=-1;

   virtual void Conj_Grad(std::vector<double>& b, std::vector<double>& x);
   virtual void jacobi(std::vector<double>& x, const std::vector<double>& b);
   virtual void Gauss_Seidel(std::vector<double>& x, const std::vector<double>& b);
   virtual void sor(std::vector<double>& x, const std::vector<double>& b, double w);


   CSRMatrix* Cholesky();
   void LU(CSRMatrix<T>& L,std::vector<T>& b_vec, std::vector<T>& x_vec);
// Private variables - there is no need for other classes
// to know about these variables
private:

};




#endif /* CSRMATRIX_H_ */
