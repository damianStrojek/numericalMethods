// Damian Strojek s184407 MN Project nr 2
#include <iostream>
#include <math.h>
#include <chrono>
#include "Matrix.h"
#define N 907
#define MAXNORM 1e-9

void residual(Matrix* A, double* b, double* x, double*& res);
double normRes(double* x);
void jacobiMethod(Matrix* A, double* b, double* x, double& norm, 
    int& iterations, double& duration);
void gaussSeidelMethod(Matrix* A, double* b, double* x, double& norm,
    int& iterations, double& duration);

int main() {
    // Task A
    int a1 = 9, a2 = -1, a3 = a2, iterations;
    Matrix* A = new Matrix(N, a1, a2, a3);
    double* b = new double[N], * x = new double[N];
    double normResidual, duration;
    for (int i = 0; i < N; i++) b[i] = sin(i * 5.0);

    // Task B
    jacobiMethod(A, b, x, normResidual, iterations, duration);
    std::cout << "\nTask B:\nJacobi method:\n\tNumber of iterations: ";
    std::cout << iterations << "\n\tDuration in milliseconds: ";
    std::cout << duration << "\n\tNorm of residual vector: ";
    std::cout << normResidual << "\n\n";

    gaussSeidelMethod(A, b, x, normResidual, iterations, duration);
    std::cout << "Gauss-Seidel method:\n\tNumber of iterations: ";
    std::cout << iterations << "\n\tDuration in milliseconds: ";
    std::cout << duration << "\n\tNorm of residual vector: ";
    std::cout << normResidual << "\n\n";

    // Garbage collector
    delete A;
    return 0;
};

// Vector residuum where res[k] = A * x[k] - b
void residual(Matrix* A, double* b, double* x, double*& res) {
    // Ready function to multiply the matrix
    res = (*A) * x;
    for (int k = 0; k < N; k++) res[k] -= b[k];
};

// Euclidean norm, formula taken from:
// https://en.wikipedia.org/wiki/Norm_(mathematics)#Euclidean_norm
double normRes(double* res) {
    double norm = 0;
    for (int i = 0; i < N; i++) norm += pow(res[i], 2);
    return sqrt(norm);
};

// Jacobi method is an iterative algorithm for determining the solutions 
// of a strictly diagonally dominant system of linear equations
void jacobiMethod(Matrix* A, double* b, double* x, double& norm, 
    int& iterations, double& duration) {
    // Setup
    iterations = 0;
    double S;
    double* res = new double[N], * xPrev = new double[N];
    for (int i = 0; i < N; i++) xPrev[i] = 1;
    
    auto start = std::chrono::steady_clock::now();
    // Main formula taken from:
    // https://en.wikipedia.org/wiki/Jacobi_method#Description
    do {
        for (int i = 0; i < N; i++) {
            S = 0;
            for (int j = 0; j < N; j++)
                if (i != j) S += A->getMatrixValueAt(i, j) * xPrev[j];
            x[i] = (b[i] - S) / A->getMatrixValueAt(i, i);
        }
        // Shifting of vector x
        for (int i = 0; i < N; i++) xPrev[i] = x[i];
        iterations++;

        residual(A, b, x, res);
        norm = normRes(res);
    } while (norm > MAXNORM);

    auto end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
};

// Gauss–Seidel method is similar to the Jacobi method. Though it can be applied 
// to any matrix with non-zero elements on the diagonals
void gaussSeidelMethod(Matrix* A, double* b, double* x, double& norm,
    int& iterations, double& duration) {
    iterations = 0;
    double S;
    double* res = new double[N], * xPrev = new double[N];
    for (int i = 0; i < N; i++) xPrev[i] = 1;

    auto start = std::chrono::steady_clock::now();
    // Main formula taken from:
    // https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method#Description
    do {
        for (int i = 0; i < N; i++) {
            S = 0;
            for (int j = 0; j < i; j++) S += A->getMatrixValueAt(i, j) * x[j];
            for (int j = i + 1; j < N; j++) S += A->getMatrixValueAt(i, j) * xPrev[j];

            x[i] = (b[i] - S) / A->getMatrixValueAt(i, i);
        }

        for (int i = 0; i < N; i++) xPrev[i] = x[i];
        iterations++;

        residual(A, b, x, res);
        norm = normRes(res);
    } while (norm > MAXNORM);

    auto end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
};