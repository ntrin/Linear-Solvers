/*
 * CSRMatrix.cpp
 *
 *  Created on: Jan 30, 2020
 *      Author: Nicolas
 */

#include <iostream>
#include "CSRMatrix.h"
#include <algorithm>
#include "math.h"
#include <vector>
#include <utility>
#include<bits/stdc++.h>
using namespace std;
#pragma once
// Constructor - using an initialisation list here
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, bool preallocate): Matrix<T>(rows, cols, false), nnzs(nnzs)
{
   // If we don't pass false in the initialisation list base constructor, it would allocate values to be of size
   // rows * cols in our base matrix class
   // So then we need to set it to the real value we had passed in
   this->preallocated = preallocate;

   // If we want to handle memory ourselves
   if (this->preallocated)
   {
      // Must remember to delete this in the destructor
      this->values = new T[this->nnzs];
      this->row_position = new int[this->rows + 1];
      this->col_index = new int[this->nnzs];
   }
}

// Constructor - now just setting the value of our T pointer
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, T *values_ptr, int *row_position, int *col_index)
: Matrix<T>(rows, cols, values_ptr), nnzs(nnzs), row_position(row_position), col_index(col_index)
{}

// destructor
template <class T>
CSRMatrix<T>::~CSRMatrix()
{
   // Delete the values array
   if (this->preallocated){
      delete[] this->row_position;
      delete[] this->col_index;
   }
   // The super destructor is called after we finish here
   // This will delete this->values if preallocated is true
}

// Explicitly print out the values in values array as if they are a matrix
template <class T>
void CSRMatrix<T>::printMatrix()
{
   std::cout << "Printing matrix" << std::endl;
   std::cout << "Values: ";
   for (int j = 0; j< this->nnzs; j++)
   {
      std::cout << this->values[j] << " ";
   }
   std::cout << std::endl;
   std::cout << "row_position: ";
   for (int j = 0; j< this->rows+1; j++)
   {
      std::cout << this->row_position[j] << " ";
   }
   std::cout << std::endl;
   std::cout << "col_index: ";
   for (int j = 0; j< this->nnzs; j++)
   {
      std::cout << this->col_index[j] << " ";
   }
   std::cout << std::endl;
}
template <class T>
void CSRMatrix<T>::printVector(std::vector<double>& v) {
	std::cout << "[ ";
	for (int i = 0; i < this->rows; i++)
		std::cout << v[i] << "  ";
	std::cout << "]" << std::endl;
}

template <class T>
CSRMatrix<T>* CSRMatrix<T>::transpose()
{
  CSRMatrix<double> *A_T = new CSRMatrix<double>(this->cols,this->rows,this->nnzs,true);
  A_T->row_position[0] = 0;

  std::vector<int> lengths(this->cols,0);
  std::vector<int> col_indices_transpose;


  for (int i = 0; i < this->rows; i++)
  {
    for (int val_index = this->row_position[i]; val_index < this->row_position[i+1]; val_index++)
    {
      int col_ind_transpose = i;
      lengths[this->col_index[val_index]]++;
      col_indices_transpose.push_back(col_ind_transpose);
    }
  }

  // for every nonzero element to be added to the transpose
//  for (int i = 0; i < A_T->rows; i++)
//  {
//    A_T->row_position[i+1] = A_T->row_position[i] + lengths[i];
//    // go through all rows in the original matrix
//    for (int val_index = this->row_position[i]; val_index < this->row_position[i+1]; val_index++)
//    {
//      A_T->values[val_index_transpose] = this->values[val_index];
//      A_T->col_index[val_index_transpose] = i;
//      val_index_transpose++;
//    }
//  }
  int val_index_transpose = 0;
  int row_find = 0;
  for (int i = 0; i < A_T->rows; i++)
  {
    A_T->row_position[i+1] = A_T->row_position[i] + lengths[i];
    // go through all nonzeros in the original matrix
    for (int val_index_2 = 0; val_index_2 < this->nnzs; val_index_2++)
    {
      if (this->col_index[val_index_2] == i)
      {
        A_T->values[val_index_transpose] = this->values[val_index_2];
        // find the row of this element
        int success = 0;
        int row_pos = 0;
        while (success == 0)
        {
          if (val_index_2 >= this->row_position[row_pos] and val_index_2 < this->row_position[row_pos+1])
          {
            A_T->col_index[val_index_transpose] = row_pos;
            success = 1;
          }
          row_pos++;
        }

//        for (int row_pos = 0; row_pos < this->rows; row_pos++)
//        {
//          if (val_index_2 >= this->row_position[row_pos])
//          {
//            A_T->col_index[val_index_transpose]
//            break;
//          }
//        }
        val_index_transpose++;
      }
    }
  }
  return A_T;
}


template <class T>
void CSRMatrix<T>::fill_random_sparse() // sparsity is between 0 and 1
{
//  srand (time(NULL));
  // initializing the arrays
  this->row_position[0] = 0; // row position array always starts at 0
  for (int i = 0; i < nnzs; i++)
  {
    this->row_position[i+1] = 0;
    this->col_index[i] = 0;
    this->values[i] = 0;
  }

  int num_diag_vals = this->rows;
  int num_non_diag_vals = this->nnzs - num_diag_vals;
  int num_non_diag_vals_per_row = num_non_diag_vals / this->rows;
  int num_vals_per_row = (num_diag_vals + num_non_diag_vals)/this->rows;
  for (int i = 0; i < this->rows; i++)
  {
    this->row_position[i+1] = (this->row_position[i] + 1);
    this->col_index[i*num_vals_per_row] = i; // column = row for diagonal
    this->values[i*num_vals_per_row] = rand() % 10 + (10*this->cols); // ensures strict diagonal dominance

    int vals_counter = 1;
    std::vector<int> vals_added(num_non_diag_vals_per_row,-1);
    while (vals_counter < num_vals_per_row)
    {
      int random_column_index = rand() % this->cols;
      if (random_column_index != i and random_column_index != (i+1)) // diagonal already filled
      {
        if (find(vals_added.begin(), vals_added.end(), random_column_index) == vals_added.end())
        {
          this->col_index[i*num_vals_per_row + vals_counter] = random_column_index;
          this->values[i*num_vals_per_row + vals_counter] = rand() % 10 + 1;
          this->row_position[i+1] += 1;
          vals_added.push_back(random_column_index);
          vals_counter++;
        }
      }
    }
  }
}


// Do a matrix-vector product
// output = this * input
template<class T>
void CSRMatrix<T>::matVecMult(std::vector<double>& input, std::vector<double>& output)
{
   if (input.size() == 0 || output.size() == 0)
   {
      std::cerr << "Input or output haven't been created" << std::endl;
      return;
   }

   // Set the output to zero
   for (int i = 0; i < this->rows; i++)
   {
      output[i] = 0.0;
   }

   int val_counter = 0;
   // Loop over each row
   for (int i = 0; i < this->rows; i++)
   {
      // Loop over all the entries in this col
      for (int val_index = this->row_position[i]; val_index < this->row_position[i+1]; val_index++)
      {
         // This is an example of indirect addressing
         // Can make it harder for the compiler to vectorise!
         output[i] += this->values[val_index] * input[this->col_index[val_index]];

      }
   }
}


// Do matrix matrix multiplication
// output = this * mat_right


// Do matrix matrix multiplication
// output = this * mat_right
template <class T>
CSRMatrix<T>* CSRMatrix<T>::matMatMult(CSRMatrix<T>& B)
{
	std::vector<int> c_row_positions;
  c_row_positions.push_back(0);
  std::vector<int> c_col_indices_main;
  int length = 0;

  for (int i = 0; i < this->rows; i++) // for all rows in the matrix product
  {
	  std::vector<int> c_col_indices;
    for (int val_index_A = this->row_position[i]; val_index_A < this->row_position[i+1]; val_index_A++)
    {
      int col_index_A = this->col_index[val_index_A];

      for (int val_index_B = B.row_position[col_index_A]; val_index_B < B.row_position[col_index_A+1]; val_index_B++)
      {
        if (find(c_col_indices.begin(), c_col_indices.end(), B.col_index[val_index_B]) == c_col_indices.end()) // to prevent double counting
        {
          c_col_indices.push_back(B.col_index[val_index_B]);
          length++;
        }

      }
    }

    // update overall structure of C
    c_col_indices_main.insert(c_col_indices_main.end(), c_col_indices.begin(), c_col_indices.end());
    c_row_positions.push_back(length);

  }

  int rows_c = c_row_positions.size()-1;
  int cols_c = *max_element(begin(c_col_indices_main), end(c_col_indices_main));
  //  int cols_c = max_element(c_col_indices_main);
  int nnzs_c = c_row_positions.back();

  CSRMatrix<double> *C = new CSRMatrix<double>(rows_c, cols_c, nnzs_c, true);

  for (int r = 0; r < c_row_positions.size(); r++)
  {
    C->row_position[r] = c_row_positions[r];
  }

  for (int v = 0; v < nnzs_c; v++)
  {
    C->col_index[v] = c_col_indices_main[v];
  }

  std::vector<double> temp_sum(C->rows,0);

  for (int i = 0; i < C->rows; i++)
  {

    for (int val_index_A_2 = this->row_position[i]; val_index_A_2 < this->row_position[i+1]; val_index_A_2++)
    {
      int col_index_A_2 = this->col_index[val_index_A_2];

      for (int val_index_B_2 = B.row_position[col_index_A_2]; val_index_B_2 < B.row_position[col_index_A_2+1]; val_index_B_2++)
      {
        // add to the temporary sum
        temp_sum[B.col_index[val_index_B_2]] += this->values[val_index_A_2] * B.values[val_index_B_2];
      }
    }

    for (int val_index_C = C->row_position[i]; val_index_C < C->row_position[i+1]; val_index_C++)
    {
//        C->values[val_index_C] = temp_sum;
      C->values[val_index_C] = temp_sum[C->col_index[val_index_C]];
      temp_sum[C->col_index[val_index_C]] = 0;
    }
  }

  C->printMatrix();
  return C;
}


template <class T>
void CSRMatrix<T>::Gauss_Seidel(std::vector<double>& x, const std::vector<double>& b)
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
      double mi;

      // Loop over all the entries in this column
      for (int val_index = this->row_position[i]; val_index < this->row_position[i+1]; val_index++)
      {

        if (i < this->col_index[val_index]) // use old x
        {
          bi -= this->values[val_index] * oldx[col_index[val_index]];
        }

        if (i == this->col_index[val_index]) // diagonal
        {
          mi = this->values[val_index];
        }

        if (i > this->col_index[val_index]) // use new x
        {
          bi -= this->values[val_index] * x[col_index[val_index]];
        }
      }
      x[i] = bi/mi;
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
void CSRMatrix<T>::sor(std::vector<double>& x, const std::vector<double>& b, double w)
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
      double mi;

      // Loop over all the entries in this column
      for (int val_index = this->row_position[i]; val_index < this->row_position[i+1]; val_index++)
      {

        if (i < this->col_index[val_index]) // use old x
        {
          bi -= this->values[val_index] * oldx[col_index[val_index]];
        }

        if (i == this->col_index[val_index]) // diagonal
        {
          mi = this->values[val_index];
        }

        if (i > this->col_index[val_index]) // use new x
        {
          bi -= this->values[val_index] * x[col_index[val_index]];
        }
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
void CSRMatrix<T>::Conj_Grad(std::vector<double>& b, std::vector<double>& x) {
	std::vector<double> r(this->rows, 0);
	std::vector<double> r1(this->rows, 0);
	std::vector<double> mult(this->rows, 0);
	std::vector<double> p(this->rows, 0);

	double tol = 0.00000001;
	int k = 0;
	double sum = 0.0;
	double dot;
	double dot3;

	this->matVecMult(x, mult);
	for (int i = 0; i < this->rows; i++){
		r[i] = b[i] - mult[i];
		p[i] = r[i];
		sum += std::abs(r[i]);
	}

	while (sum > tol){
		this->matVecMult(p, mult);
		dot = 0.0;
		dot3 = 0.0;

		for (int i = 0; i < this->rows; i++)
			dot += r[i] * r[i];

		for (int i = 0; i < this->rows; i++)
			dot3 += p[i] * mult[i];

		double alpha = dot / dot3;

		for (int i = 0; i < this->rows; i++){
			x[i] += alpha * p[i];
			r1[i] = r[i] - alpha * mult[i];
		}

		sum = 0.0;
		for (int i = 0; i < this->rows; i++)
			sum += std::abs(r1[i]);

		double beta = 0;
		for (int i = 0; i < this->rows; i++)
			beta += (r1[i] * r1[i]) / dot;
		for (int i = 0; i < this->rows; i++){
			p[i] = r1[i] + beta * p[i];
			r[i] = r1[i];
		}
		k += 1;
	}
	std::cout << "Number of iterations: " << k;
}

struct compare{ // custom comparator for sets, sort pair.second values
	            // and then sort pair.first values whilst retaining pair.second's sorting

	bool operator () (const std::pair<int, int>& left, const std::pair<int, int>& right) const {
		return left.second<right.second || (!(right.second<left.second) && left.first<right.first);
	}
};


template <class T>
void CSRMatrix<T>::jacobi(std::vector<double>& x, const std::vector<double>& b)
// Implements the Jacobi iterative solver for a sparse matrix.
{
  int iter = 0; // Keeps track of the number of iterations performed
  bool finished = false;
  double tolerance = 0.000001; // Specifies when to determine that the method has converged to a solution
  int maxiter = 1000; // Specifies when the method should "give up" (because no convergence, likely)
  while (finished == false)
  {
    std::vector<double> oldx = x;

    // Go through each row in the sparse matrix
    for (int i = 0; i < this->rows; i++)
    {
      double bi = b[i];
      double mi;

      // Go through each element in that row of the sparse matrix
      for (int val_index = this->row_position[i]; val_index < this->row_position[i+1]; val_index++)
      {
        if (i == this->col_index[val_index])
        {
          mi = this->values[val_index]; // Get diagonal element
        }
        if (i != this->col_index[val_index]) // if not on diagonal
        {
          bi -= this->values[val_index] * oldx[col_index[val_index]]; // Jacobi uses entirely values from the previous iteration for updating
        }
      }
      x[i] = bi/mi;
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
CSRMatrix<T>* CSRMatrix<T>::Cholesky() {
	//// first find how many non-zeros are in L

	std::vector<T> L_nnzs;
	std::vector<int> L_cols;
	std::vector<int> L_rows;
	auto* Ax_coord = new int[this->nnzs];

	std::set<std::pair<int, int>, compare> pairs;


	for (int i = 0; i < sizeof(this->rows) + 1; i++){ // copy row position from A, since L has similar structure
		L_rows.push_back(this->row_position[i]);
	}
	for (int k = 0; k < this->rows + 1; k++){
		for (int val_index = this->row_position[k]; val_index < this->row_position[k + 1]; val_index++){
			Ax_coord[val_index] = k; // gives x.row coordinate of the input CRS values
		}
	}

	for (int i = 0; i < this->nnzs; i++){
		if (this->col_index[i] > Ax_coord[i]){ // ignores the upper values
			for (int sub_index = Ax_coord[i] + 1; sub_index < sizeof(this->rows)+ 1; sub_index++){
				L_rows[sub_index] -= 1; // updates L_rows when a upper triangle value is ignored
			}
			continue;
		}
		else {
			L_cols.push_back(this->col_index[i]); // initialise L cols
		}
	}

	//// upper triangle ignored and L arrays updated
	// now we must account for fill in

	std::vector<int> Lx_coord;

	for (int i = 0; i < this->rows; i++){
		for (int val_index = L_rows[i]; val_index < L_rows[i + 1]; val_index++){
			Lx_coord.push_back(i); // gives x.row coordinate of the L factor
		}
	}

	typedef std::pair <int, int> IntPair;
	std::vector<IntPair> new_pairs;   // used to store the new nnzs coordinates of L
	std::vector<IntPair> tmp_pairs;   // used to evaluate the fill-ins
	// populate these pairs as we move through the nnzs

	std::vector<int> repeater;

	int counter = Lx_coord.size();
	for (int i = 0; i < counter; i++){
		auto p = std::make_pair(L_cols[i], Lx_coord[i]);
		if (std::find(tmp_pairs.begin(), tmp_pairs.end(), p) != tmp_pairs.end()){ // seeing if i pair is in the vector of pairs
			continue;
		}

		else{
			auto it = std::find(repeater.begin(), repeater.end(), L_cols[i]);
			if (it != repeater.end()) { // found it, there is a repeat
				auto index = std::distance(repeater.begin(), it); // index of the repeated col, there must be a fill-in
				auto add_p = std::make_pair(Lx_coord[index], Lx_coord[i]);

				tmp_pairs.push_back(add_p);
				new_pairs.push_back(add_p);

				tmp_pairs[index].first = -1;
				tmp_pairs[index].second = -1;
				repeater[index] = -1;
				L_cols[index] = -1;
				L_cols.push_back(add_p.first);
				Lx_coord.push_back(add_p.second);
				counter++;
			}
			tmp_pairs.push_back(p); // adds to the vector
			new_pairs.push_back(p);
			repeater.push_back(L_cols[i]);
		}
	}

	for (int i = 0; i < new_pairs.size(); i++){
		pairs.insert(new_pairs[i]);
	}

	// convert Lx_coord to L row position vector with included fill-in non-zeros

	L_rows.clear();
	L_nnzs.clear();
	L_rows.push_back(0);
	counter = 0;
	for (int j = 0; j < this->rows + 1; j++){
		for (auto const &x : pairs){
			if (j == x.second){
				counter++;
				L_nnzs.push_back(1);
			}
		}
		L_rows.push_back(counter);

	}
	for (int i = 0; i < this->rows+1; i++)
		std::cout << L_rows[i] << std::endl;
	for (int i = 0; i < L_cols.size(); i++)
	    std::cout << L_cols[i];

	delete[] Ax_coord;
}



template <class T>
CSRMatrix<T>* CSRMatrix<T>::generate_random_spd_sparse(int rows, int cols)
// Generates a random sparse matrix, with a random sparsity, that is guaranteed to be symmetric positive definite, such that it can be solved using the conjugate gradient method. Moreover, it is guaranteed to be diagonally-dominant, such that iterative methods (e.g. Jacobi, Gauss-Seidel) converge.
{


  // First, determine the structure of the sparse matrix
  std::vector<int> row_pos_vec(rows+1,0); // Stores the entries for the CSR row position array
  for (int i = 1; i <= rows; i++)
  {
    // ensures some rows have no values in them, such that the matrix is indeed quite sparse
    int random = rand() % 10;
    if (random == 0)
    {
      int vals_per_row = rand() % i;
      row_pos_vec[i] = row_pos_vec[i-1] + vals_per_row + 1;
    }
    else
      row_pos_vec[i] = row_pos_vec[i-1] + 1;

  }

  int nnzs_diag = row_pos_vec.back();

  // Create new sparse matrix with the appropriate number of nonzeros; in particular, it should be lower diagonal
  CSRMatrix<T> *low_diag = new CSRMatrix<T>(rows, cols, nnzs_diag, true);

  for (int i = 0; i <= low_diag->rows; i++)
  {
    low_diag->row_position[i] = row_pos_vec[i];
  }

  for (int i = 0; i < low_diag->rows; i++)
  {
    // First fill the diagonal elements
    low_diag->values[low_diag->row_position[i]] = rand() % 10 + (10*low_diag->rows); // Ensures matrix is strictly diagonally-dominant
    low_diag->col_index[low_diag->row_position[i]] = i;

    // Then fill the nondiagonal elements.
    std::vector<int> col_indices; // Stores the column indices for the matrix elements

    for (int val_index = low_diag->row_position[i]+1; val_index < low_diag->row_position[i+1]; val_index++)
    {
      bool added = false;
      while (!added)
      {
        int potential_col_index = rand() % i;
        if (find(col_indices.begin(), col_indices.end(), potential_col_index) == col_indices.end()) // Ensures that values are not added to the same position in the matrix
        {
          low_diag->col_index[val_index] = potential_col_index;
          low_diag->values[val_index] = rand() % 10 + 1;
          col_indices.push_back(potential_col_index);
          added = true;
        }
      }
    }
    col_indices.clear();
  }

  // The product of a lower triangular matrix with nonzero diagonal elements, with its transpose (an upper triangular matrix, also with nonzero diagonal elements) is an SPD matrix.
  CSRMatrix<T> *upper_diag = low_diag->transpose();
  CSRMatrix<T> *spd_sparse = low_diag->matMatMult(*upper_diag); // Note that thanks to polymorphism, matMatMult here is from specifically the CSR matrix class, rather than the generic Matrix class

  return spd_sparse;
}




template <class T>
CSRMatrix<T>* CSRMatrix<T>::generate_random_diagdom_sparse(int rows, int cols, double sparsity)
// fills a sparse-formatted (old Yale format) matrix with random values. Note that the number of non-zeros in the matrix has already been specified in the constructor.
{
  srand ( time(NULL) );

  // determines number of nonzeros to add to the sparse matrix. Note that since this method generates a matrix with a uniform number of nonzero elements in each row, it adjusts the user specified sparsity accordingly.
  double nnzs_input = (rows * cols) * sparsity;
  nnzs = round(nnzs_input / static_cast<double>(rows)) * static_cast<double>(rows);

  CSRMatrix<T> *sparse_mat = new CSRMatrix<T>(rows, cols, nnzs, true);

  // initializing the arrays for row position, column index, and value
  sparse_mat->row_position[0] = 0; // row position array always starts at 0
  for (int i = 0; i < nnzs; i++)
  {
    sparse_mat->row_position[i+1] = 0;
    sparse_mat->col_index[i] = 0;
    sparse_mat->values[i] = 0;
  }
  // Sets parameters for filling the matrix elements
  int num_diag_vals = sparse_mat->rows;
  int num_non_diag_vals = sparse_mat->nnzs - num_diag_vals;
  int num_non_diag_vals_per_row = num_non_diag_vals / sparse_mat->rows;
  int num_vals_per_row = (num_diag_vals + num_non_diag_vals)/sparse_mat->rows;
//std::cout << "num vals" << num_vals_per_row;
  // Go through all rows in the sparse matrix
  for (int i = 0; i < sparse_mat->rows; i++)
  {
    // First add the diagonal entry
    sparse_mat->row_position[i+1] = (sparse_mat->row_position[i] + 1); // each row has at least one element, the diagonal entry
    sparse_mat->col_index[i*num_vals_per_row] = i; // column = row for diagonal entry
    sparse_mat->values[i*num_vals_per_row] = rand() % 10 + (10*sparse_mat->cols); // Ensures strict diagonal dominance
    int vals_counter = 1; // Stores the number of values added to the matrix row thus far
    std::vector<int> vals_added(num_non_diag_vals_per_row,-1); // Stores the entries to add to this matrix row
    // Fills the given number of nonzero elements for every row, with random indices and random values
    while (vals_counter < num_vals_per_row)
    {
      int random_column_index = rand() % sparse_mat->cols;
      if (random_column_index != i and random_column_index != (i+1)) // Diagonal is already filled
      {
        if (find(vals_added.begin(), vals_added.end(), random_column_index) == vals_added.end()) // If value has not already been added to that column position in the row, add it
        {
          sparse_mat->col_index[i*num_vals_per_row + vals_counter] = random_column_index;
          sparse_mat->values[i*num_vals_per_row + vals_counter] = rand() % 10 + 1;
          sparse_mat->row_position[i+1] += 1;
          vals_added.push_back(random_column_index);
          vals_counter++;

        }
      }
    }
  }


  return sparse_mat;
}


template <class T>
void CSRMatrix<T>::matRowVecMult(std::vector<double>& input, std::vector<double>& output)
// Multiplies a ROW vector by a matrix, outputting the row vector solution
{
	std::vector<double> temp_sums(input.size(),0); // to add together non-consecutive element multiplications

  // go through all rows in the matrix
  for (int i = 0; i < this->rows; i++)
  {
    // go through all elements in that row
    for (int val_index = this->row_position[i]; val_index < this->row_position[i+1]; val_index++)
    {
      // multiply by corresponding vector entry and add to temp sum
      temp_sums[this->col_index[val_index]] += this->values[val_index] * input[this->col_index[val_index]];
    }
  }

  // copy into output vector
  for (int i = 0; i < output.size(); i++)
  {
    output[i] = temp_sums[i];
  }

}

template <class T>
void CSRMatrix<T>::LU(CSRMatrix<T>& L, std::vector<T>& b_vec,  std::vector<T>& x_vec)
{
   std::cout << "LU SPARSE decomposition"<< std::endl;
   //declare value n
   int n = this->cols;
   //copy of the original b vector to use in the LU forward sustitutiom part
   std::vector<T> b_vec_L_copy(n);
   b_vec_L_copy = b_vec;
   //SPARSE
   // ---------------------------------------------------
   //1.FIND DIAGONAL ELEMENTs
   //initialise diagonal elements to store them
   std::vector<T> diag(n);
   for (int i = 0 ; i < n ; i++)
   {diag[i] = 0; cout << "diag[" << i << "]= " << diag[i] << endl;}
   //Iterate through every row and compare to column position to find the coincident (i.e.) Diagonal elements
   for(int i = 0; i < n ; i++)
   {
      //Compute Number of elements per row
      int n_el_row = this->row_position[i+1] - this->row_position[i];
      //Find the diagonal element
      for (int j=0; j< n_el_row ; j++)
      {
         if (this->col_index[this->row_position[i] + j ]  ==  i ){
            diag[i] = this->values[this->row_position[i]+j];
            cout << "Diagonal element diag[" << i << "]= " << diag[i] << endl;
         }
      }

      //Initialise auxiliary variables to do gaussian elimination
      double d;
      double multiplo = 0.0;
      //3.DIAGONALISATION
      //Loop over each column to eliminate nonzero entries. Here the variable kk is used to scrape every position of the entire values vector
      //to try to find the values whose column is the one in which we are trying to eliminate any nonzero entry.
      for (int kk=this->row_position[i+1] ; kk < this->nnzs; kk++)
      {
         //Get index of whichever nnz entry column value concides with column currently we are in ( i-th. column)
         if (this->col_index[kk] == i)
         {
            //Multiple compute to vanish elements in i-th column
            d = diag[i];
            //set to double
            cout << d/1.0 << endl;
            //Compute the alfa factor, here called multiplo to do row operations
            multiplo =  this->values[kk]/d;
            //While statement to find the row at which that nonzero entry at column i belongs to.
            //This is used then to access the directly to row_position vector and use the value stored there as index to then operate.
            int r_update = 0;
            while (this->row_position[r_update] < kk) // this->row_position[r_update] < kk
            {r_update +=1;}

            cout << "ROW ------ >   ------ >  " << r_update << endl;


            cout << "KK   : " << kk <<endl;
            cout << "multiplo ---> col: " << i << ", row: "<< r_update << " = " << this->values[kk]/d << endl;
            cout << "Row in: " << r_update << endl;
            //Initialisation of the variables col_pivot and col_update
            //This responds to the requeriment to know the column at which they belong the rest of elements
            //of the rows in pivot row and row being updated. These are crutial to do row operation updates ------------------------
            //---------------------------------------------> The MOST COMPLEX PART OF LU and GAUSSIAN ELIM.
            int col_pivot = 0;
            int col_update= 0;
            //4. ROW OPERATIONS : UPDATE VALUES OF ROW
            //----------------------------------------
            //Updating the values of any given operation at any row will result in either 2 cases, ramins the same, because there is zeros in the pivot row at j-col
            //or its updated. Then if its updated there would be 2 cases, that the nnzs ofrow being updated are the same ( ie cols of pito rows coincide with row_updated)
            //1. which will require of doing -----> row_updated[j] -= mult* pivot_row[j] (j-th. element in pivot row, which does not represent the column )
            //2. There is no coincidence and do not exist an entry at row_updated[j] which coincides with the nnz entry from pivot_row. --> update (( this->VECTORS  ))

            for (int r= this->row_position[i]; r < this->row_position[i+1];r++)
            {
               cout << " ----------------------------------------------" << endl;
               //counter to activate if NEW nnz appears
               int activ =0;
               //Acces the column of the nnz entries of row[i];
               col_pivot = this->col_index[r];


               //loop over the nnz entries of the update row[r_update] --> ABSOLUTELY! (NBA - Ricky rubio! #SPAIN-#Barcelona-#Masnou)
               // 2 cases (activ == 0 || 1 )
               for (int rr= this->row_position[r_update]+i; rr< this->row_position[r_update+1]; rr++)
               {
                  //Iterate with row r_pivot
                  col_update = this->col_index[rr];
                  //2 cases
                  if (col_pivot == col_update)
                  {this->values[rr] -= multiplo*this->values[r];activ +=1;}
                  else{activ +=0;}
               }
               for (int jjh = 0; jjh < this->nnzs; jjh++)
               {
                  cout << "Values updated row (no new nnz) [" << jjh << "] :: " << this->values[jjh] << endl;
               }
               // 5. INCLUDE A NEW NON ZERO ENTRY TO THE UPDATED ROW by changing all three vectors that characterise CSRMatrix
               if (activ==0)
               {
                  //Visualisation of matrix values prior to the
                  for (int h = 0; h < this->nnzs; h++)
                  {cout << "Values prior row update [" << h << "] :: " << this->values[h] << endl;}
                  //This methodology expoits the condition that we have previously allocated a larger amount of memory than actually needed by seting a nnzs greater
                  //INITIALISE nnzs and nnzs_Filled (which are the oroginal nnz values of the sparse matrix)
                  //new_nonzero entry position to fill in the nnzs array ( to keep order)
                  ////cout << "col pivot --> " << col_pivot << endl;
                  int new_nnz_r = 0;
                  for (int lim= this->row_position[r_update]; lim< this->row_position[r_update+1]; lim++)
                  {
                     //update the new_nnz_r positioner so that is ordered;
                     if (this->col_index[lim] < col_pivot){
                     new_nnz_r +=1;}
                  }
                  // 5.1 - Push back all of the entries to the right from the back and allocate new nnz entry
                  for(int nnz_up= this->nnzs; nnz_up >= (this->row_position[r_update]+new_nnz_r); nnz_up--)
                  {
                     this->values[nnz_up] = this->values[nnz_up-1];
                  }

                  for (int h = 0; h < this->cols+1; h++)
                  {
                     cout << "Sanity check row_pos [" << h << "] = " << this->row_position[h] << endl;
                  }

                  cout << "new_nnz_index --> " << new_nnz_r << endl;


                  this->values[this->row_position[r_update]+new_nnz_r] = 0.0 -  multiplo*this->values[r]; // r -- this->row_position[i] + new_nnz_r
                  for (int jh = 0; jh < this->nnzs; jh++)
                  {
                     cout << "Values [" << jh << "] :: " << this->values[jh] << endl;
                  }
                  cout << "-------------------------------------" << endl;
                  cout << "New entry --> "  << 0.0 -  multiplo*this->values[this->row_position[i] + new_nnz_r]  << endl;
                  // 5.2 - Update the columns and fix the new nnz entry column
                  for(int nnz_up=this->nnzs+1; nnz_up > (this->row_position[r_update]+new_nnz_r); nnz_up--){this->col_index[nnz_up] = this->col_index[nnz_up-1];}
                  this->col_index[this->row_position[r_update]+col_pivot] = col_pivot;
                  // 5.3 - Update n#elements in the current row_update & rest of rows(rp_up ==== row_position_updater)
                  for (int rp_up=0; rp_up < (this->rows-r_update); rp_up++){this->row_position[r_update+1 + rp_up] +=1;}
               }
               for (int cc = 0; cc < this->nnzs; cc++)
               {
                  cout << "COLUMN POSITION [" << cc << "] = " << this->col_index[cc] << endl;
               }
               cout << "------ UPDATED ROW ------------ " << r_update << "-th. row" << endl;
            }
            //STORE MULT VALUES IN THE CORRESPONDANT POSITION IN A SPARSE LOWER TRIANGULAR MATRIX
            //L.values (( L.values.size =   [[ n + 2*(n-1) ]]  ))
            int lr =0;
            while (this->row_position[lr] < kk) {lr +=1;} //get the row we are in
            L.values[L.row_position[lr]] = multiplo; //access and set values to the multiple used to perform element elimination in gauss elim.
            cout << "L values being printed " << endl;
            for (int i=0 ; i < L.col_index[n+1]; i++) {cout << L.values[i] << endl;}
         }

      }
      cout << "" << endl;
      cout << "FULL CYCLE COLUMNA  -----> " << i+1 << "_a  columna  ----------------------------------  F U L L CYCLE    -   C O L U M N : " << i <<  "  !   " << endl;
      cout << " - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl;
   }
   cout << "---- all finished -------------------------------------------------------------------------------" << endl;
}




























