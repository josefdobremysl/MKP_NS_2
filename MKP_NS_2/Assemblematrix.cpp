#include "Assemblematrix.h"

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <armadillo>

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <tuple>
#include <cmath>
#include <numeric>

using namespace Eigen;
using namespace std;
using namespace arma;

SparseMatrixVectorResult AssembledMatrix(vector<int>& I, vector<int>& J, vector<double>& VAL, vector<double>& bD) {
	SparseMatrixVectorResult assembledMatrix;

	// Initialize sparse matrix and right-hand side vector

	Eigen::SparseMatrix<double> A(bD.size(), bD.size());
	// Assemble matrix and right-hand side vector directly
	for (size_t k = 0; k < I.size(); ++k) {
		A.coeffRef(I[k], J[k]) += VAL[k]; // Accumulate values for duplicate indices
	}

	A.makeCompressed();

	// Assign assembled matrix and right-hand side vector to the result structure
	assembledMatrix.sparsematrix = A;
	assembledMatrix.right_side_vector = bD;



	return assembledMatrix;
}