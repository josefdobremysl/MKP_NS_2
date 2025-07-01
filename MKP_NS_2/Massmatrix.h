#ifndef MASSMATRIX_H
#define MASSMATRIX_H

#include "mesh_processor.h"

#include <string>
#include <vector>
#include <tuple>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

SparseMatrixVectorResult MassMatrix(const string& meshFile, MatrixVectorResult& meshResult);

#endif // MASSMATRIX_H
