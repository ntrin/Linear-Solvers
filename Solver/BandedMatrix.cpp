/*
 * BandedMatrix.cpp
 *
 *  Created on: Feb 2, 2020
 *      Author: Nicolas
 */
#include "BandedMatrix.h"
#include <iostream>
#include <vector>
#include "math.h"
using namespace std;


// Constructor - using an initialisation list here
// ml and mu represents the number of elements beyond a tridiagonal matrix in the lower and upper triangular parts, respectively.
template <class T>
BandedMatrix<T>::BandedMatrix(int rows, int cols, int ml, int mu, bool preallocate): Matrix<T>(rows, cols, false), ml(ml), mu(mu), bw(ml + mu + 1), num_vals((bw+2) * rows)
{
  this->preallocated = preallocate;
  // If we want to handle memory ourselves
  if (this->preallocated)
  {
    this->values = new T[this->num_vals];
  }
}
template <class T>
BandedMatrix<T>::BandedMatrix(int rows, int cols, int ml, int mu, T *values_ptr): Matrix<T>(rows, cols, values_ptr)
{}
template <class T>
BandedMatrix<T>::~BandedMatrix()
{}
template <class T>
void BandedMatrix<T>::printMatrix()
{
  cout << endl << "Banded matrix values (including placeholder zeros at the top left and bottom right): " << endl;
  for (int i = 0; i < this->rows; i++)
  {
    for (int j = 0; j < this->bw+2; j++)
      cout << this->values[i*(this->bw + 2) + j] << " ";
    cout << endl;
  }
  cout << endl;
}
template<class T>
void BandedMatrix<T>::matVecMult(vector<double>& input, vector<double>& output)
// Multiplies a banded matrix by a column vector, storing the result in a vector output
{
  // Loop over each row in the matrix
  for (int i = 0; i < this->rows; i++)
  {

    int j = 0; // Stores the column position in banded format
    for (int col_index = i - (this->ml + 1); col_index <= i + (this->mu + 1); col_index++) // col_index gives starting index for column vector multiplication
    {

      if (col_index >= 0 and col_index < this->cols) // Note that col_index < 0 means the multiplication will be excluded since the matrix values are placeholders
      {
        output[i] += this->values[i * (this->bw + 2) + j] * input[col_index];
      }
      j++;
    }
  }
}

template <class T>
void BandedMatrix<T>::fill_random_banded()
{

  for (int i = 0; i < this->rows; i++)
  {
    int j = 0;
    for (int col_index = i - (this->ml + 1); col_index <= i + (this->mu + 1); col_index++)
    {
      if (col_index >= 0 and col_index < this->cols)
      {
        if (col_index == i) // if element is on the diagonal
          this->values[i * (this->bw + 2) + j] = rand() % 10 + (10*this->cols);
        if (col_index != i) // if element is not on diagonal
//          else
        {
          this->values[i * (this->bw + 2) + j] = rand() % 10 + 1;
        }
      }
      j++;
    }
  }
}


template <class T>
void BandedMatrix<T>::jacobi(vector<double>& x, const vector<double>& b)
// Implement the Jacobi solver for a banded matrix
{
  int iter = 0; // Keeps track of the number of iterations performed
  double tolerance = 1e-8; // Specifies when to determine that the method has converged to a solution
  int maxiter = 1000; // Specifies when the method should "give up" (because no convergence, likely)
  bool finished = false;
  while (finished == false)
  {
    vector<double> oldx = x;
    // Loop over each row
    for (int i = 0; i < this->rows; i++)
    {
      double bi = b[i];
      double aii;
      int j = 0;
      // Loop over elements in bandwidth on that row
      for (int col_index = i - (this->ml + 1); col_index <= i + (this->mu + 1); col_index++)
      {
        if (col_index >= 0 and col_index < this->cols)
        {
          if (col_index == i) // if element is on the diagonal
            aii = this->values[i * (this->bw + 2) + j];
          if (col_index != i) // if element is not on diagonal
//          else
          {
            bi -= this->values[i * (this->bw + 2) + j] * oldx[col_index];
          }
        }
        j++;
      }
      x[i] = bi/aii;
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
//    cout << "norm difference = " << norm_diff << endl;
    if (iter++ == maxiter)
    {
      finished = true;
      cout << "Failed to converge after " << maxiter << " iterations..." << endl;
    }
    if (abs(norm_diff) < tolerance)
    {
      finished = true;
      cout << "Iterations for convergence = " << iter << endl;
    }
  }
}

template <class T>
void BandedMatrix<T>::Gauss_Seidel(std::vector<double>& x, const std::vector<double>& b) //don't want to modify m and b, but do want to modify x.
// implement for a sparse matrix
{
  int iter = 0;
  int success = 0;
  double tolerance = 0.000001;
  int maxiter = 1000;
  while (success == 0)
  {
    std::vector<double> oldx = x;
    // Loop over each row
    for (int i = 0; i < this->rows; i++)
    {
      double bi = b[i];
      double aii;
      int j = 0;

      // Loop over elements in bandwidth on that row
      for (int col_index = i - (this->ml + 1); col_index <= i + (this->mu + 1); col_index++)
      {

        if (col_index >= 0 and col_index < this->cols)
        {
          if (col_index == i) // if element is on the diagonal
            aii = this->values[i * (this->bw + 2) + j];
          if (col_index < i) // if element is not on diagonal
          {
            bi -= this->values[i * (this->bw + 2) + j] * x[col_index];
          }
          if (col_index > i) // if element is not on diagonal
          {
            bi -= this->values[i * (this->bw + 2) + j] * oldx[col_index];
          }
        }
        j++;
      }
      x[i] = bi/aii;
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
//    cout << "norm difference = " << norm_diff << endl;
    if (iter++ == maxiter)
    {
      success = 1;
      std::cout << "Failed to converge after " << maxiter << " iterations..." << std::endl;
    }
    if (abs(norm_diff) < tolerance)
    {
      success = 1;
      std::cout << "Iterations for convergence = " << iter << std::endl;
    }
  }
}

template <class T>
void BandedMatrix<T>::sor(std::vector<double>& x, const std::vector<double>& b, double w) //don't want to modify m and b, but do want to modify x.
// implement for a sparse matrix
{
  int iter = 0;
  int success = 0;
  double tolerance = 1e-7;
  int maxiter = 1000;
  while (success == 0)
  {
    std::vector<double> oldx = x;
    // Loop over each row
    for (int i = 0; i < this->rows; i++)
    {
      double bi = b[i];
      double aii;
      int j = 0;

      // Loop over elements in bandwidth on that row
      for (int col_index = i - (this->ml + 1); col_index <= i + (this->mu + 1); col_index++)
      {

        if (col_index >= 0 and col_index < this->cols)
        {
          if (col_index == i) // if element is on the diagonal
            aii = this->values[i * (this->bw + 2) + j];
          if (col_index < i) // if element is not on diagonal
//          else
          {
            bi -= this->values[i * (this->bw + 2) + j] * x[col_index];
          }
          if (col_index > i) // if element is not on diagonal
//          else
          {
            bi -= this->values[i * (this->bw + 2) + j] * oldx[col_index];
          }
        }
        j++;
      }
      x[i] = bi/aii;
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
    // to test if converged
    double norm_diff = norm_old - norm_new;
    if (iter++ == maxiter)
    {
      success = 1;
      std::cout << "Failed to converge after " << maxiter << " iterations..." << std::endl;
    }
    if (abs(norm_diff) < tolerance)
    {
      success = 1;
      std::cout << "Iterations for convergence = " << iter << std::endl;
    }
  }
}






