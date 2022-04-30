// Damian Strojek s184407 MN Project nr 2
#pragma once

class Matrix {
private:
	int size;
	double** A;
public:
	Matrix(int N, int a1, int a2, int a3);
	Matrix(const Matrix& M);
	~Matrix();
	double** setMatrix(double** A) { this->A = A; };
	int setSize(const int size) { this->size = size; };
	void setAAt(const int i, const int j, double value) {
		this->A[i][j] = value;
	};
	double** getMatrix() { return this->A; };
	// I wanted to make `size` and `A` private and the operator[] didn't cooperate
	double getAAt(const int i, const int j) { 
		return this->A[i][j]; 
	};
	int getSize() { return this->size; };
	void writeMatrix();
	Matrix& operator=(const Matrix &M);
	double* operator*(const double *v);
};
