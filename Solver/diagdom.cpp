/*
 * main_solver.cpp
 *
 *  Created on: Jan 23, 2020
 *      Author: Nicolas
 */
/*
#include "Matrix.cpp"
#include "Matrix.h"
#include "Gauss_Elim.cpp"
#include "Gauss_Elim.h"

void makediagdom(){
	int rows = 10;
	int cols = 10;
	bool wantdiagdom = false;
	bool sparse = false;
	int min = 0;
	int max = 10;
	// to generate random integers within a defined range
	std::random_device rd;     // only used once to initialise (seed) engine
	std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
	std::uniform_int_distribution<int> uni(min,max);

	auto random_integer = uni(rng);


	auto *matrix_data = new double [rows * cols];
//	auto *mat = new Matrix(rows, cols, true);
//	auto *mat = new Matrix<double>(rows, cols, matrix_data);

	for (int i = 0; i < rows * cols; i++) // fill matrix with random integers
		matrix_data[i] = random_integer;




	std::cout << "This program generates a square matrix and completes any (?) kind of solver method.\n";
	std::cout << "How many rows should the square matrix have? ";
	std::cin >> rows;
	cols = rows;

	std::cout << "Would you like a diagonal dominant matrix? true/false: ";
	std::cin >> wantdiagdom;
	if (wantdiagdom == true){
		for (int i = 0; i < rows; i++){ // iterate through rows
			int sum_row = 0; // sum of rows
			for (int j = 0; j < cols; j++) // iterate through cols
				sum_row += std::abs(matrix_data[i + j]);

			while (matrix_data[i * cols + i] <= std::abs(sum_row - matrix_data[i * cols + i])){ // check if diagonal entry is the largest in the row
				matrix_data[i * cols + i] += 1; // add to the entry if not true
			}
		}

	}

	std::cout << "Would you like a sparse matrix? true/false: ";
	std::cin >> sparse; // then we would need to use the CSRMatrix class.

	delete[] matrix_data;
}

*/



