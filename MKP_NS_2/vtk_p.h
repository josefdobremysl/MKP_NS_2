#ifndef VTK_P_H
#define VTK_P_H

#include "mesh_processor.h"

#include <string>
#include <vector>
#include <tuple>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

void VTKfile_p(const string& filename, const string& meshFile, MatrixVectorResult& meshResult, Eigen::VectorXd& solution_x);

#endif // VTK_P_H