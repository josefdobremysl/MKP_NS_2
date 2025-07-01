#ifndef NEUMANN_H
#define NEUMANN_H

#include <string>
#include <vector>
#include <tuple>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;


vector<double> NeumannBoundaryCondition(const string& meshFile, MatrixVectorResult& meshResult, const int& tag, const double& t, vector<double> w);

#endif // NEUMANN_H