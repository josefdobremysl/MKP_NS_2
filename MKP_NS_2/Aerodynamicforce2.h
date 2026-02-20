#ifndef AERODYNAMICFORCE2_H
#define AERODYNAMICFORCE2_H

#include <string>
#include <vector>
#include "mesh_processor.h"

// Vrací {Fx, Fy} z integrálu trakce po hranách s daným tagem.
std::vector<double> AerodynamicForce2(
    const std::string& meshFile,
    MatrixVectorResult& meshResult,
    const int& tag,
    const double& t,
    const std::vector<double>& vx_vy_p
);

#endif
