#ifndef VTK_H
#define VTK_H

#include "mesh_processor.h"

#include <string>
#include <vector>
#include <tuple>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

void VTKfile(const string& filename, const string& meshFile, MatrixVectorResult& meshResult, Eigen::VectorXd& solution_x);

#endif // DIRICHLET_H
