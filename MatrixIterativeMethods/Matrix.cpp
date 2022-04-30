// Damian Strojek s184407 MN Project nr 2
#pragma once
#include <iostream>
#include "Matrix.h"

// Creating new matrix NxN where main diagonal values = a1, 
// surrounding values of diagonals = a2, and then a3
Matrix::Matrix(int N, int a1, int a2, int a3) {
	this->A = new double* [N];
	this->size = N;

	// Creating new matrix using two dimensional array
	for (int i = 0; i < N; i++) this->A[i] = new double[N];
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == j) this->A[i][j] = a1;
			else if (i - 1 == j || i + 1 == j) this->A[i][j] = a2;
			else if (i - 2 == j || i + 2 == j) this->A[i][j] = a3;
			else this->A[i][j] = 0;
		}
	}
};

// Copy of Matrix M
Matrix::Matrix(const Matrix& M) {
	this->size = M.size;
	this->A = new double* [M.size];

	for (int i = 0; i < this->size; i++) this->A[i] = new double[this->size];
	for (int i = 0; i < this->size; i++)
		for (int j = 0; j < this->size; j++) this->A[i][j] = M.A[i][j];
};

// Ability to debug the matrix (pretty big write)
void Matrix::writeMatrix() {
	std::system("cls");
	for (int i = 0; i < this->size; i++) {
		std::cout << i << " [ ";
		for (int j = 0; j < this->size; j++)
			// To look good for integers, didnt do it for floats
			if (this->A[i][j] >= 0) std::cout << " " << this->A[i][j] << " ";
			else std::cout << this->A[i][j] << " ";
		std::cout << "]\n";
	}
};

Matrix& Matrix::operator=(const Matrix& M) {
	for (int i = 0; i < this->size; i++)
		for (int j = 0; j < this->size; j++) this->A[i][j] = M.A[i][j];
	return *this;
};

double* Matrix::operator*(const double* v) {
	double* u = new double[this->size];

	for (int i = 0; i < this->size; i++) {
		double temp = 0;
		for (int j = 0; j < this->size; j++) temp += this->A[i][j] * v[j];
		u[i] = temp;
	}
	return u;
};

Matrix::~Matrix() {
	for (int i = 0; i < this->size; i++) delete this->A[i];
	delete this->A;
}
