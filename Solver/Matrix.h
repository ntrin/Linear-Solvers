/*
 * Matrix.h
 *
 *  Created on: Jan 23, 2020
 *      Author: Nicolas
 */

#ifndef MATRIX_H_
#define MATRIX_H_
#pragma once
#include <vector>

template <class T>
class Matrix{
public:
	Matrix(int rows, int cols, bool preallocate); // allocate ourselves
    Matrix(int rows, int cols, T *values_ptr); // allocated from outside

    virtual ~Matrix();

    void printValues();
	virtual void printMatrix();
	virtual void printVector(std::vector<double>& v);

    T *values = nullptr;
    int rows = -1;
    int cols = -1;

    virtual void matVecMult(std::vector<double>& input, std::vector<double>& output);
    virtual void matRowVecMult(std::vector<double>& input, std::vector<double>& output);
    void Transposer(Matrix<T>& transposed);


    void Gauss_Elim(std::vector<double>& b, std::vector<double>& x);
    void backsubstitution(std::vector<double>& b, std::vector<double>& x);
    void Forward_Substitution(std::vector<double>& b_vec , std::vector<double>& y_vec);
    void FS(Matrix<T>& Lower, std::vector<T>& y_vec , std::vector<T>& b_vec); // this version accepts matrix as argument
    void BS(std::vector<T>& x_vec, std::vector<T>& y_vec);

    // Note that while all three matrix classes have a Jacobi method,
    // due to polymorphism this function is specific to the dense Matrix class
    virtual void jacobi(std::vector<T>& x, std::vector<T>& b);
    virtual void Gauss_Seidel(std::vector<T>& x, std::vector<T>& b);
    virtual void sor(std::vector<double>& x, const std::vector<double>& b, const double w);


    virtual void Conj_Grad(std::vector<T>& b, std::vector<T>& x);

    // DIRECT METHODS
    void Cholesky(Matrix<T>&L, std::vector<T>& b, std::vector<T>& x); // returns pointer to L, cholesky factor
    void LU(Matrix<T>& L,std::vector<T>& b_vec, std::vector<T>& x_vec);

    void fill_random_dense();

    Matrix<T>* generate_random_spd_dense(int rows, int cols);

    void matMatMult(Matrix<T>& mat_right, Matrix<T>& output);

    void fill_spd_dense();
    void fill_random_dense_diagdom(); // fills dense matrix with random values


protected:
    bool preallocated = false;

private:
    int size_of_values = -1;
};




#endif /* MATRIX_H_ */
