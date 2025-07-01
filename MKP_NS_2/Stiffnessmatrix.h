#ifndef STIFFNESSMATRIX_H
#define STIFFNESSMATRIX_H

#include "mesh_processor.h"

#include <string>
#include <vector>
#include <tuple>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

SparseMatrixVectorResult StiffnessMatrix(const string& meshFile, MatrixVectorResult& meshResult, Eigen::VectorXd x, const double& t);

#endif // STIFFNESSMATRIX_H
