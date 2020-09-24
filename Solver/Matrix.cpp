/*
 * Matrix.cpp
 *
 *  Created on: Jan 23, 2020
 *      Author: Nicolas
 */
#include <iostream>
#include "math.h"
#include <memory>
#include <random>
#include "Matrix.h"
#include <time.h>
#pragma once

template <class T>
Matrix<T>::Matrix(int rows, int cols, bool preallocate): rows(rows), cols(cols), size_of_values(rows * cols), preallocated(preallocate)
{
   // If we want to handle memory ourselves
   if (this->preallocated)
   {
      // Must remember to delete this in the destructor
      this->values = new T[size_of_values];
   }
}


template <class T>
Matrix<T>::Matrix(int rows, int cols, T *values_ptr)
	: rows(rows), cols(cols), size_of_values(rows * cols), values(values_ptr){};

// destructor
template <class T>
Matrix<T>::~Matrix(){
   // Delete the values array
   if (this->preallocated){
      delete[] this->values;
   }
};

// Just print out the values in our values array
template <class T>
void Matrix<T>::printValues() {
   std::cout << "Printing values" << std::endl;
	for (int i = 0; i< this->size_of_values; i++) {
      std::cout << this->values[i] << " ";
   }
   std::cout << std::endl;
}

// Explicitly print out the values in values array as if they are a matrix
template <class T>
void Matrix<T>::printMatrix() {
   std::cout << "Printing matrix" << std::endl;
   for (int j = 0; j < this->rows; j++){
      std::cout << std::endl;
      for (int i = 0; i < this->cols; i++){
         // We have explicitly used a row-major ordering here
         std::cout << this->values[i + j * this->cols] << " ";
      }
   }
   std::cout << std::endl;
}

template <class T>
void Matrix<T>::printVector(std::vector<double>& v) {
	std::cout << "[ ";
	for (int i = 0; i < this->rows; i++)
		std::cout << v[i] << "  ";
	std::cout << "]" << std::endl;
}

template <class T>
void Matrix<T>::matMatMult(Matrix<T>& mat_right, Matrix<T>& output)
{
  // Check our dimensions match
  if (this->cols != mat_right.rows)
  {
   std::cerr << "Input dimensions for matrices don't match" << std::endl;
   return;
  }
  // Check if our output matrix has had space allocated to it
  if (output.values != nullptr)
  {
   // Check our dimensions match
   if (this->rows != output.rows || this->cols != output.cols)
   {
     std::cerr << "Input dimensions for matrices don't match" << std::endl;
     return;
   }
  }
  // The output hasn't been preallocated, so we are going to do that
  else
  {
   output.values = new T[this->rows * mat_right.cols];
  }
  // Set values to zero before hand
  for (int i = 0; i < output.size_of_values; i++)
  {
   output.values[i] = 0;
  }
  // Now we can do our matrix-matrix multiplication
  // CHANGE THIS FOR LOOP ORDERING AROUND
  // AND CHECK THE TIME SPENT
  // Does the ordering matter for performance. Why??
  for(int i = 0; i < this->rows; i++)
  {
   for(int k = 0; k < this->cols; k++)
   {
     for(int j = 0; j < mat_right.cols; j++)
     {
        output.values[i * output.cols + j] += this->values[i * this->cols + k] * mat_right.values[k * mat_right.cols + j];
     }
   }
  }
}

template <class T>
void Matrix<T>::matVecMult(std::vector<double>& input, std::vector<double>& output)
{
  for (int i = 0; i < this->rows; i++)
  {
    output[i] = 0; // ensures output vector initialized to zero.
    for (int j = 0; j < this->cols; j++)
    {
      output[i] += this->values[i*this->rows + j] * input[j];
    }
  }
}

template <class T>
void Matrix<T>::matRowVecMult(std::vector<double>& input, std::vector<double>& output)
{
  for (int i = 0; i < this->cols; i++)
  {
    output[i] = 0; // ensures output vector initialized to zero.
    for (int j = 0; j < this->rows; j++)
    {
      output[i] += this->values[j*this->cols + i] * input[j];
    }
  }
}

template <class T>
void Matrix<T>::Transposer(Matrix<T>& transposed) {
  for (int i = 0; i < this->rows; i++){
    for (int j = 0; j < this->cols; j++){ // <= so that diagonals are nonzero
      transposed.values[i*this->cols + j] = this->values[j*this->cols + i];
    }
  }
}

template <class T>
void Matrix<T>::Gauss_Elim(std::vector<double>& b, std::vector<double>& x) {
	for (int i = 0; i < this->rows - 1; i++){
		double diag = this->values[i * this->cols + i]; // diagonal

		for (int j = i + 1; j < this->rows; j++){
			double lambda = this->values[j * this->cols + i] / diag;
			// multiply whole i row elements by lambda and subtract to j row cols element
			for (int k = 0; k < this->rows; k++)
				this->values[j*this->cols + k] -= lambda * this->values[i*this->cols + k];
			b[j] -= lambda * b[i];
		}
	}
}

template <class T>
void Matrix<T>::backsubstitution(std::vector<double>& b, std::vector<double>& x) {
	for (int i = this->rows - 1; i >= 0; i--){
		double tmp = b[i];
		for (int j = this->rows - 1; j > i; j--)
			tmp -= x[j] * this->values[i*this->cols + j];
		x[i] = tmp / this->values[i*this->cols + i];
	}
}

template <class T>
void Matrix<T>::Forward_Substitution(std::vector<double>& b_vec , std::vector<double>& y_vec)
{
   //From the Upper matrix and the vector b we can obtain the solution
   // 1. FORWARD_SUBSTITUTION
   ////  Solve type problems --> L(Ux) = b ; L(y) = b x////
   // ¡¡ REMEMBER TO SET INITIAL VALUE TO THE VECTOR b in MAIN to DO WELL !! // b_vec = {8,6,4,2};
   for (int i = 0; i <this->rows; i++){
      double tmp = b_vec[i];
      for (int j = 0 ; j < i; j++)
         tmp -= y_vec[j] * this->values[i*this->rows+ j];
      y_vec[i] = tmp/ this->values[i*this->rows + i];
   }
}

template <class T>
void Matrix<T>::jacobi(std::vector<T>& x, std::vector<T>& b)
// Implements the Jacobi iterative solver for a dense matrix.
{
  int iter = 0; // Keeps track of the number of iterations performed
  int n = this->rows;
  double tolerance = 0.000001; // Specifies when to determine that the method has converged to a solution
  int maxiter = 1000; // Specifies when the method should "give up" (because no convergence, likely)

  bool finished = false;
  while (finished == false)
  {
	  std::vector<T> oldx = x;

    // Go through each row in the dense matrix
    for (int i = 0; i < n; i++)
    {
      double bi = b[i];
      double aii = this->values[i * cols + i]; // Get diagonal entry for that row

      // Go through each element in that row of the dense matrix
      for (int j = 0; j < n; j++)
      {
        if (j != i)
          bi -= this->values[j + cols * i] * oldx[j]; // Jacobi uses entirely values from the previous iteration for updating
      }
      x[i] = bi/aii; // Update the x vector value corresponding to that row
    }

    // Compute norm of the previous iteration's x vector
    double total_old = 0;
    for (int i = 0; i < oldx.size(); i++)
    {
      total_old += oldx[i] * oldx[i];
    }
    double norm_old = sqrt(total_old);

    // Compute norm of the new iteration's x vector
    double total_new = 0;
    for (int i = 0; i < x.size(); i++)
    {
      total_new += x[i] * x[i];
    }
    double norm_new = sqrt(total_new);

    // Compare the norms to determine whether the method has converged to a solution
    double norm_diff = norm_old - norm_new;

    if (iter++ == maxiter)
    {
      finished = true;
      std::cout << "Failed to converge after " << maxiter << " iterations..." << std::endl;
    }

    if (abs(norm_diff) < tolerance)
    {
      finished = true;
      std::cout << "Iterations for convergence = " << iter << std::endl;
    }
  }
}


template <class T>
void Matrix<T>::Gauss_Seidel(std::vector<T>& x, std::vector<T>& b)
// Implements the Gauss-Seidel iterative solver for a dense matrix.
{
  int iter = 0; // Keeps track of the number of iterations performed
  int n = this->rows;
  double tolerance = 1e-12; // Specifies when to determine that the method has converged to a solution
  int maxiter = 10000; // Specifies when the method should "give up" (because no convergence, likely)

  bool finished = false;
  while (finished == false)
  {
    std::vector<T> oldx = x;

    // Go through each row in the dense matrix
    for (int i = 0; i < n; i++)
    {
      double bi = b[i];
      double aii = this->values[i * cols + i]; // Get diagonal entry for that row

      for (int j = 0; j < i; j++)
      {
        bi -= this->values[j + cols * i] * x[j]; // For the lower triangular portion of the matrix, the Gauss-Seidel method uses the already computed, new iteration's x vector
      }

      for (int j = i+1; j < n; j++)
      {
        bi -= this->values[j + cols * i] * oldx[j]; // For the lower triangular portion of the matrix, the Gauss-Seidel method uses the previous iteration's x vector
      }
      x[i] = bi/aii;
    }

    // Compute norm of the old iteration's x vector
    double total_old = 0;
    for (int i = 0; i < oldx.size(); i++)
    {
      total_old += oldx[i] * oldx[i];
    }
    double norm_old = sqrt(total_old);

    // Compute norm of the new iteration's x vector
    double total_new = 0;
    for (int i = 0; i < x.size(); i++)
    {
      total_new += x[i] * x[i];
    }
    double norm_new = sqrt(total_new);

    // Compare the norms to determine whether the method has converged to a solution
    double norm_diff = norm_old - norm_new;
    if (iter++ == maxiter)
    {
      finished = true;
      std::cout << "Failed to converge after " << maxiter << " iterations..." << std::endl;
    }

    if (abs(norm_diff) < tolerance)
    {
      finished = true;
      std::cout << "Iterations for convergence = " << iter << std::endl;
    }
  }
}


template <class T>
void Matrix<T>::sor(std::vector<double>& x, const std::vector<double>& b, const double w)
// Implements the Successive Over-Relaxation iterative solver for a dense matrix.
{
  int iter = 0;
  int n = rows;
  double tolerance = 0.000001;
  int maxiter = 1000;
  bool finished = false;
  while (finished == false)
  {
	  std::vector<double> oldx = x;
    for (int i = 0; i < n; i++) // for each row
    {
      double bi = b[i];
      double mi = this->values[i * cols + i];

      for (int j = 0; j < i; j++)
      {
        bi -= this->values[j + cols * i] * x[j];
      }

      for (int j = i+1; j < n; j++)
      {
        bi -= this->values[j + cols * i] * oldx[j];
      }
      x[i] = bi/mi;
      x[i] *= w; // acceleration parameter; w=1 recovers Gauss-Seidel
      x[i] += (1-w)*oldx[i];
    }

    // compute norm of old vector
    double total_old = 0;
    for (int i = 0; i < oldx.size(); i++)
    {
      total_old += oldx[i] * oldx[i];
    }
    double norm_old = sqrt(total_old);

    // compute norm of new vector
    double total_new = 0;
    for (int i = 0; i < x.size(); i++)
    {
      total_new += x[i] * x[i];
    }
    double norm_new = sqrt(total_new);

    double norm_diff = norm_old - norm_new;

    if (iter++ == maxiter)
    {
      finished = true;
      std::cout << "Failed to converge after " << maxiter << " iterations..." << std::endl;
    }

    if (abs(norm_diff) < tolerance)
    {
      finished = true;
      std::cout << "Iterations for convergence = " << iter << std::endl;
    }
  }
}

template <class T>
void Matrix<T>::Conj_Grad(std::vector<T>& b, std::vector<T>& x) {

	std::vector<double> r(this->rows, 0); // initialize some vectors necessary for the calculation
	std::vector<double> r1(this->rows, 0);
	std::vector<double> mult(this->rows, 0);
	std::vector<double> p(this->rows, 0);

	double tol = 1e-1; // set tolerance
	int k = 0; // initialize number of iterations
	double sum = 0.0; // initialize variables
	double dot;
	double dot3;

	this->matVecMult(x, mult);
	for (int i = 0; i < this->rows; i++){ // initialize r
		r[i] = b[i] - mult[i];
		p[i] = r[i];
		sum += std::abs(r[i]);
	}

	while (sum > tol){
		this->matVecMult(p, mult);
		dot = 0.0;
		dot3 = 0.0;

		for (int i = 0; i < this->rows; i++) // first, calculate the scalars needed for the alpha coefficient
			dot += r[i] * r[i];

		for (int i = 0; i < this->rows; i++)
			dot3 += p[i] * mult[i];

		double alpha = dot / dot3; // alpha coefficient

		for (int i = 0; i < this->rows; i++){ // update x and r for the possible next iterations
			x[i] += alpha * p[i];
			r1[i] = r[i] - alpha * mult[i];
		}

		sum = 0.0;
		for (int i = 0; i < this->rows; i++) // if the r values are small enough to satisfy our tolerance, exit the loop
			sum += std::abs(r1[i]);

		double beta = 0;
		for (int i = 0; i < this->rows; i++) // beta coefficient
			beta += (r1[i] * r1[i]) / dot;
		for (int i = 0; i < this->rows; i++){
			p[i] = r1[i] + beta * p[i]; // update p and r values
			r[i] = r1[i];
		}
		k++;

	}

	std::cout << "Number of iterations: " << k << std::endl;

}

template <class T>
void Matrix<T>::Cholesky(Matrix<T>& L, std::vector<T>& b, std::vector<T>& x) {
//    auto* L = new T[this->rows * this->cols];

	for (int i = 0; i < this->rows; i++) // cycle through the rows of A
		for (int j = 0; j < i + 1; j++){ // but only the lower triangle
			double s = 0;
			for (int k = 0; k < j; k++) // iterate to find s values
				s += L.values[i * this->rows + k] * L.values[j * this->rows + k];
			if (i == j) // at diagonal s becomes square of L [i,j]
				L.values[i * this->rows + j] = sqrt(this->values[i * this->rows + i] - s);
			else
				L.values[i * this->rows + j] = 1.0 / L.values[j * this->rows + j] * (this->values[i * this->rows + j] - s);
		}
	for (int i = 0; i < this->rows * this->cols; i++)
		if (L.values[i] > 1e10 || (L.values[i] < 0.0001 && L.values[i] > 0)) // compiler gives values on the order of 1e-317 due to rounding errors which can be approximated to 0.
			L.values[i] = 0;


	std::vector<double> y (this->rows, 0); // initialize vectors needed for back, forward substitution
    double* LT_values = new double[this->rows*this->cols]; // initialize empty L transpose matrix
    for (int i = 0; i < this->rows * this->cols; i++)
    	LT_values[i] = 0;
                                                                 // initialize matrices for the Cholesky factor
    auto* LT = new Matrix<T>(this->rows, this->cols, LT_values); // and the corresponding transpose

    L.Transposer(*LT);
    L.Forward_Substitution(b, y);
    LT->backsubstitution(y, x); // The result is x
}


template <class T>
void Matrix<T>::fill_random_dense()
{
  for (int i = 0; i < this->rows * this->cols; i++)
  {
    if (i % (this->rows+1) == 0)
      this->values[i] = rand() % 10 + (9*this->rows); // ensures strictly diagonally dominant
    else
      this->values[i] = rand() % 10;
  }
}



template <class T>
Matrix<T>* Matrix<T>::generate_random_spd_dense(int rows, int cols)
// Generates a random dense matrix, of dimensions specified in the argument,
{
  std::srand ( time(NULL) );

  // Create new sparse matrix with the appropriate number of nonzeros; in particular, it should be lower diagonal
  Matrix<double> *A = new Matrix<double>(rows, cols, true);
  A->fill_random_dense();

  // The product of a lower triangular matrix with nonzero diagonal elements, with its transpose (an upper triangular matrix, also with nonzero diagonal elements) is an SPD matrix.
  Matrix<double> *B = new Matrix<double>(rows, cols, true);
  B->fill_random_dense();

  Matrix<double> *spd_dense = new Matrix<double>(rows, cols, true);
  A->matMatMult(*B, *spd_dense);

  return spd_dense;
}




template <class T>
void Matrix<T>::fill_spd_dense()
{
  // first create a random lower triangular matrix, with diagonal entries non-zero
  for (int i = 0; i < this->rows; i++)
  {
    for (int j = 0; j < this->cols; j++) // <= so that diagonals are nonzero
    {
      if (j <= i)
      {
        this->values[i*this->cols + j] = rand() % 10 + 1;
      }
      else
        this->values[i*this->cols + j] = 0;
    }
  }

//  // then create its transpose
  Matrix<double> *transpose = new Matrix<double>(this->rows, this->cols, true);
  for (int i = 0; i < this->rows; i++)
  {
    for (int j = 0; j < this->cols; j++) // <= so that diagonals are nonzero
    {
      transpose->values[i*this->cols + j] = this->values[j*this->cols + i];
    }
  }


  Matrix<double> *output = new Matrix<double>(this->rows, this->cols, true);

  // the product of these two matrices will be SPD
  this->matMatMult(*transpose, *output);

  for (int i = 0; i < this->rows; i++)
  {
    for (int j = 0; j < this->cols; j++) // <= so that diagonals are nonzero
    {
      this->values[i*this->cols + j] = output->values[i*this->cols + j];
    }
  }
}

template <class T>
void Matrix<T>::fill_random_dense_diagdom()
{
  for (int i = 0; i < this->rows * this->cols; i++)
  {
    if (i % (this->rows+1) == 0)
      this->values[i] = rand() % 10 + (9*this->rows); // ensures strictly diagonally dominant
    else
      this->values[i] = rand() % 10;
  }
}


template <class T>
void Matrix<T>::FS(Matrix<T>& Lower, std::vector<T>& y_vec , std::vector<T>& b_vec)
{
   //From the Upper matrix and the vector b we can obtain the solution
   // 1. FORWARD_SUBSTITUTION
   ////  Solve type problems --> L(Ux) = b ; L(y) = b ////
   // ¡¡ REMEMBER TO SET INITIAL VALUE TO THE VECTOR b in MAIN to DO WELL !! // b_vec = {8,6,4,2};
   for (int i = 0; i <rows; i++){
      double tmp = b_vec[i];
      for (int j = 0 ; j < i; j++)
         tmp -= y_vec[j] * Lower.values[i*rows+ j];
        y_vec[i] = tmp / this->values[i*cols + i];
      std::cout << "y[" << i << "]= " << y_vec[i]<< std::endl;
    }
}



template <class T>
void Matrix<T>::LU(Matrix<T>& L, std::vector<T>& b_vec,  std::vector<T>& x_vec)
{
    std::cout << "LU decomposition"<< std::endl;
   //declare value n
   int n = this->cols;
   //copy of the original b vector to use in the LU forward sustitutiom part
   std::vector<T> b_vec_L_copy(n);
   b_vec_L_copy = b_vec ;
   //Gauss Elimination - Find an upper triang matrix from a given
    for (int i = 0; i < n - 1; i++)
    {
      // diagonal elements
        double diag = this->values[i*this->cols + i];
        L.values[i*this->cols + i] = 1;
        for (int j = i + 1; j < n; j++)
        {
            double mult = this->values[i + j * this->cols] / diag; // iterate through matrix elements going downwards
            this->values[i + j * this->cols] -= mult * diag;
            L.values[i + j * this->cols] = mult;
            //Update elements in the row
            for (int k=i+1; k<n; k++){
            this->values[k + j*this->cols] -= mult* this->values[k + i*this->cols];}
            //Update the constant vector b
            b_vec[j] -= mult * b_vec[i];
        }
    }
    L.values[n*n-1]= 1;
   // TWO-STAGE SOLVE Ax = b  --> LU x = b --> L(Ux) = b ;
   // -----------------------------------------------------------------------------------------------------
   // 1. Ly = b -- FS
   // ---------------
   //declare vector y to store values
   std::vector<T> vector_y(n);
   FS(L, vector_y , b_vec_L_copy);
   // 2. Ux = y -- BS
   // ---------------
   BS(x_vec, vector_y);
   // -----------------------------------------------------------------------------------------------------
}


template <class T>
void Matrix<T>::BS(std::vector<T>& x_vec, std::vector<T>& y_vec)
{
   //From the Upper matrix and the vector b we can obtain the solution
   ////  solve type problems --> Ux = y ////
   for (int i = rows - 1; i >= 0; i--){
      double tmp2 = y_vec[i];
      for (int j = rows - 1; j > i; j--)
         tmp2 -= x_vec[j] * this->values[i*cols + j];
        x_vec[i] = tmp2 / this->values[i*cols + i];
    }

}
































