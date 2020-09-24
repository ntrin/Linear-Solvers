/*
 * BandedMatrix.h
 *
 *  Created on: Feb 2, 2020
 *      Author: Nicolas
 */

#ifndef BANDEDMATRIX_H_
#define BANDEDMATRIX_H_
#pragma once
#include "Matrix.h"

template <class T>
class BandedMatrix: public Matrix<T>
{
public:
  // constructor where we want to preallocate ourselves; ml and mu indicate the number of elements to the left of the tridiagonal band, respectively
  BandedMatrix(int rows, int cols, int ml, int mu, bool preallocate);

  // constructor where we already have allocated memory outside
  BandedMatrix(int rows, int cols, int ml, int mu, T *values_ptr);

  // destructor
  ~BandedMatrix();
  // Print out the values in our matrix
  virtual void printMatrix();
  // Perform some operations with our matrix
//  void matMatMult(DiagonalMatrix<T>& mat_right, DiagonalMatrix<T>& output);
  void matVecMult(std::vector<double>& input, std::vector<double>& output);
  void fill_random_banded();

  void jacobi(std::vector<double>& x, const std::vector<double>& b); // Solves banded matrix with Jacobi method
  void Gauss_Seidel(std::vector<double>& x, const std::vector<double>& b); // Solves banded matrix with Gauss-Seidel method
  void sor(std::vector<double>& x, const std::vector<double>& b, double w); // Solves banded matrix with Successive Over-Relaxation method


  int ml = 0;
  int mu = 0;
  int bw = 0; // bandwidth
  int num_vals = 0;


// Private variables - there is no need for other classes
// to know about these variables
private:

};




#endif /* BANDEDMATRIX_H_ */
