#ifndef AERODYNAMICFORCE_H
#define AERODYNAMICFORCE_H

#include <string>
#include <vector>
#include <tuple>
#include <Eigen/Sparse>
#include "mesh_processor.h"

using namespace Eigen;
using namespace std;


vector<double> AerodynamicForce(const string& meshFile, MatrixVectorResult& meshResult, const int& tag, const double& t, vector<double> vx_vy_p);

#endif // AERODYNAMICFORCE_H
