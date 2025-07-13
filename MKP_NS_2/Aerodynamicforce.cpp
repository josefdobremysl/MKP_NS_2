#include "mesh_processor.h"
#include "Dirichlet.h"
#include "Neumann.h"
#include "Aerodynamicforce.h"
#include "function.h"

#include <Eigen/Sparse>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <tuple>
#include <cmath>
#include <armadillo>
#include <numeric>
#include <omp.h>

using namespace Eigen;
using namespace std;
using namespace arma;

vector<double> AerodynamicForce(const string& meshFile, MatrixVectorResult& meshResult, const int& tag, const double& t, vector<double> vx_vy_p) {

	double rho = 1.25;
	double F_x = 0;
	double F_y = 0;

	MatrixVectorResult loadedMesh = meshResult;
	int nTri = loadedMesh.integerResult;		// poèet trojúhelníkù
	int nNod = loadedMesh.matrix3.size();
	vector<int> M(nNod, 0);
	for (int k = 0; k < nTri; k++) {
		int A = loadedMesh.matrix1[k][0] - 1;					// Oznaèení vrcholù trojúhelníku k
		int B = loadedMesh.matrix1[k][1] - 1;
		int C = loadedMesh.matrix1[k][2] - 1;
		vector<int> verTri = { A, B, C };
		for (int i = 0; i < 3; i++) {
			M[verTri[i]] = 1;
		}
	}
	int nVert = 0;

	for (int k = 0; k < nNod; k++) {
		nVert += M[k];

	}

	vector<int> vectorOfTags = loadedMesh.vector1;
	vector<vector<int>> lineVertex = loadedMesh.matrix2;
	int size = loadedMesh.integerResult2;
	vector<vector<int>> markedEdges(size, vector<int>(2, 0));

	vector<vector<double>> nM(size, vector<double>(2, 0));
	int nNodes = 2 * loadedMesh.matrix3.size() + nVert;
	int N = loadedMesh.matrix3.size();
	vector<double> F_x_inPoints(loadedMesh.matrix3.size(), 0);
	vector<double> F_y_inPoints(loadedMesh.matrix3.size(), 0);

	for (int i = 0; i < size; i++) {
		if (vectorOfTags[i] == tag) {
			int D = lineVertex[i][0] - 1;
			int E = lineVertex[i][1] - 1;
			int F = lineVertex[i][2] - 1;

			vector<double> cooD = loadedMesh.matrix3[D];		// Souøadnice vrcholù trojúhelníku k
			vector<double> cooE = loadedMesh.matrix3[E];

			double pD = vx_vy_p[D + 2 * N];
			double pE = vx_vy_p[E + 2 * N];

			vector<double> nn = { cooE[1] - cooD[1], cooD[0] - cooE[0] };
			double n_norm = (sqrt(nn[0] * nn[0] + nn[1] * nn[1]));
			vector<double> n = { nn[0] / n_norm,nn[1] / n_norm };

			//////////////////////////////////////////////////////////////////////////////

			double Q1 = 0.5 - sqrt(3.0 / 5.0) / 2.0;
			double Q2 = 0.5;
			double Q3 = 0.5 + sqrt(3.0 / 5.0) / 2.0;


			vector<double> cooF1 = { cooD[0] + (cooE[0] - cooD[0]) * Q1, cooD[1] + (cooE[1] - cooD[1]) * Q1 };
			vector<double> cooF2 = { cooD[0] + (cooE[0] - cooD[0]) * Q2, cooD[1] + (cooE[1] - cooD[1]) * Q2 };
			vector<double> cooF3 = { cooD[0] + (cooE[0] - cooD[0]) * Q3, cooD[1] + (cooE[1] - cooD[1]) * Q3 };

			double pF1 = pD + (pE - pD) * Q1;
			double pF2 = pD + (pE - pD) * Q2;
			double pF3 = pD + (pE - pD) * Q3;

			double edgeLength = sqrt(pow(cooD[0] - cooE[0], 2) + pow(cooD[1] - cooE[1], 2));


			vector<double> p_F = { pF1, pF2, pF3 };

			double w = 1;

			vector<double> nx = { n[0],n[0], n[0] };
			vector<double> ny = { n[1],n[1], n[1] };

			F_x_inPoints[D] += rho * edgeLength * w * W3dotprod(nx, p_F);
			F_x_inPoints[E] += rho * edgeLength * w * W3dotprod(nx, p_F);
			F_x_inPoints[F] += rho * edgeLength * w * W3dotprod(nx, p_F);
			F_y_inPoints[D] += rho * edgeLength * w * W3dotprod(ny, p_F);
			F_y_inPoints[E] += rho * edgeLength * w * W3dotprod(ny, p_F);
			F_y_inPoints[F] += rho * edgeLength * w * W3dotprod(ny, p_F);


		}
	}


	for (int i = 0; i < F_x_inPoints.size(); i++) {
		F_x += F_x_inPoints[i];
		F_y += F_y_inPoints[i];

	}


	vector<double> F = { -F_x,-F_y };

	return F;

};