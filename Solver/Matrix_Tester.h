/*
 * Matrix_Tester.h
 *
 *  Created on: Feb 2, 2020
 *      Author: Nicolas
 */

#ifndef MATRIX_TESTER_H_
#define MATRIX_TESTER_H_
#pragma once
//#include "Matrix.h"
#include "Matrix.cpp"
//#include "CSRMatrix.h"
#include "CSRMatrix.cpp"

template <class T>
class Matrix_Tester {
public:

	void Gauss_Elim_test(int rows, int cols, std::vector<T>& b_val);
	void jacobi_test(int rows, int cols, bool sparsed, bool banded, std::vector<T>& b_val);
    void Gauss_Seidel_test(int rows, int cols, bool sparsed, bool banded, std::vector<T>& b_val);
    void SOR_test(int rows, int cols, bool sparsed, bool banded, std::vector<T>& b_val);
	void CG_test(int rows, int cols, bool sparsed, std::vector<T>& b_val);
	void Cholesky_test(int rows, int cols, std::vector<T>& b_val);
	bool SPD_dense_test(Matrix<T>& candidate_spd_dense_mat);

	bool sparsed = false;
	bool banded = false;
private:

};




#endif /* MATRIX_TESTER_H_ */
