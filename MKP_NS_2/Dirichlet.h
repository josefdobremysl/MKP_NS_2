#pragma once
#ifndef DIRICHLET_H
#define DIRICHLET_H

#include "mesh_processor.h"

#include <string>
#include <vector>
#include <tuple>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

SparseMatrixVectorResult DirichletBoundaryCondition(const string& meshFile, MatrixVectorResult& meshResult, const int& tag, const double& t, vector<double> w);

#endif // DIRICHLET_H