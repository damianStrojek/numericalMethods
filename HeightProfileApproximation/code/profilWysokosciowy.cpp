#include <iostream>
#include <math.h>
#include <string>
#include <fstream>
#include <chrono>
#include "Matrix.h"
// Reference samples - each .txt has more than 500 lines
#define SAMPLES 500

// Samples are defined as distance of the route from point 0
// and elevation above sea level
struct sample{
	double x;
	double y;
};

void pivoting(Matrix& U, Matrix& L, Matrix& P, int i);
void lowerUpperDecomposition(Matrix* A, double* b, double* x, int N);
double polynomialLagrange(const sample* samples, double distance, const int n);
double splines(const sample* samples, double distance, const int nodes_number);
void interpolation(const sample* samples, const int METHOD, const int nodes_number,
	sample* nodes, const char* filename, double* duration);
bool read_data(const char* path, sample* samples);

int main(){
	const char* data_sets[4] = { "genoa_rapallo.txt", "ostrowa.txt", "diff_heights.txt", "tczew_starogard.txt" };
	sample* samples = new sample[SAMPLES];	
	int intervals[4] = { 16, 40, 60, 80 };
	const int METHOD = 1;
	double average_duration[4] = { 0 };

	for (int i = 0; i < 4; ++i){
		const char* filename = data_sets[i];
		std::string directory_path = "C:\\Users\\Administrator\\source\\repos\\profilWysokosciowy\\data\\";
		std::string file_path = directory_path;
		file_path.append(filename);
		filename = file_path.c_str();

		if (!read_data(filename, samples)) return -1;

		for (int j = 0; j < 4; ++j){
			int step = intervals[j];

			const int nodes_number = SAMPLES / step;
			sample* nodes = new sample[nodes_number];

			for (int k = 0, l = 0; l < nodes_number; k += step, l++) {
				nodes[l].x = samples[k].x;
				nodes[l].y = samples[k].y;
			}

			// method = 0 LAGRANGE, method = 1 SPLINES
			interpolation(samples, METHOD, nodes_number, nodes, filename, &average_duration[j]);
		}
	}

	if (!METHOD) std::cout << "\n Average duration(Lagrange method) :\n";
	else if(METHOD) std::cout << "\n Average duration(splines method) :\n";

	for (int i = 0; i < 4; ++i){
		int nodes_number = SAMPLES / intervals[i];
		double duration = average_duration[i] / 4;
		std::cout << "\t" << nodes_number <<" nodes: " << duration << " ms\n";
	}

	return 0;
};

void pivoting(Matrix& U, Matrix& L, Matrix& P, int i){
	// Finding appropriate pivot
	double pivot = abs(U.getAAt(i, i));
	int pivot_idx= i;

	for (int j = i + 1; j < U.getSize(); j++)
		if (abs(U.getAAt(j, i)) > pivot){
			pivot = abs(U.getAAt(j, i));
			pivot_idx = j;
		}

	if (!U.getAAt(pivot_idx, i)){
		std::cout << "ERROR: Matrix is singular.\n";
		return;
	}

	// Interchanging rows
	if (pivot_idx != i){
		double tmp;

		for (int j = 0; j < U.getSize(); j++) {
			if (j >= i) {
				tmp = U.getAAt(i, j);
                U.setAAt(i, j, U.getAAt(pivot_idx, j));
                U.setAAt(pivot_idx, j, tmp);
			}
			else {
				tmp = L.getAAt(i, j);
                L.setAAt(i, j, L.getAAt(pivot_idx, j));
                L.setAAt(pivot_idx, j, tmp);
			}

            tmp = P.getAAt(i, j);
            P.setAAt(i, j, P.getAAt(pivot_idx, j));
            P.setAAt(pivot_idx, j, tmp);
		}
	}
};

// LU Decomposition factors a matrix as the product of a 
// lower triangular matrix and an upper triangular matrix
void lowerUpperDecomposition(Matrix* A, double* b, double* x,
    int N){
    // Upper and lower triangular matrixes
    Matrix* L = new Matrix(N), * U = new Matrix(*A), * P = new Matrix(N);

    for(int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            if (i == j) {
                L->setAAt(i, j, 1);
                P->setAAt(i, j, 1);
            }

    // Create L and U, such that A = L * U
    for (int i = 0; i < N - 1; i++){
        pivoting(*U, *L, *P, i);

        for (int j = i + 1; j < N; j++){
            L->setAAt(j, i,
                (U->getAAt(j, i) / U->getAAt(i, i)));
            for (int k = i; k < N; k++)
                U->setAAt(j, k,
                    (U->getAAt(j, k) - L->getAAt(j, i) * U->getAAt(i, k)));
        }
    }

    // Solving L * y = b for y going forward
    b = (*P) * b;
    double* y = new double[N];
    double S;
    for (int i = 0; i < N; i++){
        S = 0;
        for (int j = 0; j < i; j++) S += L->getAAt(i, j) * y[j];
        y[i] = (b[i] - S) / L->getAAt(i, i);
    }

    // Solving U * x = y going backward
    for (int i = N - 1; i >= 0; i--){
        S = 0;
        for (int j = i + 1; j < N; j++) S += U->getAAt(i, j) * x[j];
        x[i] = (y[i] - S) / U->getAAt(i, i);
    }

    delete L; delete U; delete P;
};

// For each point we calculate polynomial
double polynomialLagrange(const sample* samples, const double distance, const int n){
	double temp = 0;

    for (int i = 0; i < n; ++i){
        double a = 1;

        for (int j = 0; j < n; ++j)
            if (i != j) a *= (distance - samples[j].x) / (samples[i].x - samples[j].x);

        temp += a * samples[i].y;
    }

    return temp;
};

//
double splines(const sample* samples, double distance, const int nodes_number) {
	int N = 4 * (nodes_number - 1); // (n-1) equations
	Matrix* M = new Matrix(N);
	double* b = new double[N], * x = new double[N];

	for (int i = 0; i < N; i++) {
		b[i] = 0;
		x[i] = 1;
	}

	// Building system of equations
	M->setAAt(0, 0, 1);
	b[0] = samples[0].y;

	double h;

	h = samples[1].x - samples[0].x;
	M->setAAt(1, 0, 1);
	M->setAAt(1, 1, h);
	M->setAAt(1, 2, pow(h, 2));
	M->setAAt(1, 3, pow(h, 3));
	b[1] = samples[1].y;

	M->setAAt(2, 2, 1);
	b[2] = 0;

	h = samples[nodes_number - 1].x - samples[nodes_number - 2].x;
	M->setAAt(3, 4 * (nodes_number - 2) + 2, 2);
	M->setAAt(3, 4 * (nodes_number - 2) + 3, 6 * h);
	b[3] = 0;

	for (int i = 1; i < nodes_number - 1; i++) {
		h = samples[i].x - samples[i - 1].x;

		M->setAAt(4 * i, 4 * i, 1);
		b[4 * i] = samples[i].y;

		M->setAAt(4 * i + 1, 4 * i, 1);
		M->setAAt(4 * i + 1, 4 * i + 1, h);
		M->setAAt(4 * i + 1, 4 * i + 2, pow(h, 2));
		M->setAAt(4 * i + 1, 4 * i + 3, pow(h, 3));
		b[4 * i + 1] = samples[i + 1].y;

		M->setAAt(4 * i + 2, 4 * (i - 1) + 1, 1);
		M->setAAt(4 * i + 2, 4 * (i - 1) + 2, 2 * h);
		M->setAAt(4 * i + 2, 4 * (i - 1) + 3, 3 * pow(h, 2));
		M->setAAt(4 * i + 2, 4 * i + 1, -1);
		b[4 * i + 2] = 0;

		M->setAAt(4 * i + 3, 4 * (i - 1) + 2, 2);
		M->setAAt(4 * i + 3, 4 * (i - 1) + 3, 6 * h);
		M->setAAt(4 * i + 3, 4 * i + 2, -2);
		b[4 * i + 3] = 0;
	}

	lowerUpperDecomposition(M, b, x, N);

	double elevation{};
	for (int i = 0; i < nodes_number - 1; i++) {
		elevation = 0;

		if (distance >= samples[i].x && distance <= samples[i + 1].x) {
			for (int j = 0; j < 4; j++) {
				double h = distance - samples[i].x;

				elevation += x[4 * i + j] * pow(h, j);
			}
			break;
		}
	}

	delete[] x; delete[] b; delete M;
	return elevation;
};

// mode = 0 LAGRANGE, mode = 1 SPLINES
void interpolation(const sample* samples, const int METHOD, const int nodes_number,
	sample* nodes, const char* path, double* duration) {
	// setting up the file
	std::string file_path = path;

	if (!METHOD) file_path.replace(55, 4, "results\\Lagrange");
	else if (METHOD) file_path.replace(55, 4, "results\\splines");
	file_path.replace(file_path.end() - 4, file_path.end(), "_");
	file_path.append(std::to_string(nodes_number));
	file_path.append(".txt");

	std::ofstream file;
	file.open(file_path, file.out);

	if (!file) {
		std::cout << "Failed to save the results (" << nodes_number << ").\n";
		return;
	}

	double result{};
	auto start = std::chrono::high_resolution_clock::now();

	for (double i = nodes[0].x; i <= nodes[nodes_number - 1].x; i += 8) {
		bool interpolate = true;

		for (int j = 0; j < nodes_number; j++) {
			// values at given nodes don't need to be interpolated
			if ((int)nodes[j].x == i) {
				interpolate = false;
				result = nodes[j].y;
				break;
			}
		}

		if (interpolate) {
			if (!METHOD) result = polynomialLagrange(nodes, i, nodes_number);
			else if (METHOD) result = splines(nodes, i, nodes_number);
		}

		file << i << " " << result << "\n";
	}

	file.close();

	auto end = std::chrono::high_resolution_clock::now();
	auto difference = end - start;
	*duration += std::chrono::duration<double, std::milli>(difference).count();

	file_path.replace(file_path.end() - 4, file_path.end(), "_nodes.txt");

	file.open(file_path, file.out);

	if (!file.is_open()) {
		std::cout << "Failed to save the results of nodes.txt (" << nodes_number << ").\n";
		return;
	}

	for (int i = 0; i < nodes_number; i++) file << nodes[i].x << " " << nodes[i].y << "\n";

	file.close();
};

// Read data from source files and save them in samples
bool read_data(const char* filename, sample* samples)
{
	std::ifstream file;
	file.open(filename, file.in);

	if (!file.is_open()){
		std::cout << "Failed to open the specified file: " << filename << "\n";
		return false;
	}

	std::string sample;

	for (int i = 0; i < SAMPLES; i++){
		getline(file, sample, ' ');
		samples[i].x = std::stod(sample);
		getline(file, sample);
		samples[i].y = std::stod(sample);
	}

	file.close();
	return true;
};