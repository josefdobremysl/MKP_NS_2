#ifndef ASSEMBLEMATRIX_H
#define ASSEMBLEMATRIX_H

#include "mesh_processor.h"

#include <string>
#include <vector>

#include <Eigen/Sparse>


using namespace Eigen;
using namespace std;

SparseMatrixVectorResult AssembledMatrix(vector<int>& I, vector<int>& J, vector<double>& VAL, vector<double>& bD);

#endif // ASSEMBLEMATRIX_H