// Damian Strojek s184407 MN Projekt nr 3
#pragma once
#include <iostream>
#include <math.h>
#include "Matrix.h"

// Creating new matrix NxN
Matrix::Matrix(int N) {
	this->A = new double* [N];
	this->size = N;

	// Creating new matrix using two dimensional array
	for (int i = 0; i < N; i++) this->A[i] = new double[N];
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			this->A[i][j] = 0;
};

Matrix::Matrix(const Matrix& M)
{
	this->size = M.size;
	this->A = new double* [M.size];

	for (int i = 0; i < this->size; i++)
		this->A[i] = new double[this->size];

	for (int i = 0; i < this->size; i++)
		for (int j = 0; j < this->size; j++)
			this->A[i][j] = M.A[i][j];
};

// Ability to debug the matrix (pretty big write)
void Matrix::writeMatrix() {
	std::system("cls");
	for (int i = 0; i < this->size; i++) {
		std::cout << i << " [ ";
		for (int j = 0; j < this->size; j++)
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
	for (int i = 0; i < this->size; i++) delete[] this->A[i];
	delete[] this->A;
};