#include "mesh_processor.h"
#include "Dirichlet.h"
#include "Neumann.h"
#include "function.h"
#include "Assemblematrix.h"

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <armadillo>
#include <omp.h>

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <tuple>
#include <cmath>
#include <numeric>

using namespace Eigen;
using namespace std;
using namespace arma;


SparseMatrixVectorResult DirichletBoundaryCondition(const string& meshFile, MatrixVectorResult& meshResult, const int& tag, const double& t, vector<double> w) {


	SparseMatrixVectorResult resultS;
	MatrixVectorResult result1 = meshResult;
	vector<int> vectorOfTags = result1.vector1;
	int nTri = result1.integerResult;		// po�et troj�heln�k�D
	int nNodes = result1.matrix3.size();			// po�et bod�
	int nNodes2 = 2 * nNodes * 16;
	vector<int> I(nNodes2, 0);
	vector<int> J(nNodes2, 0);
	vector<double> VAL(nNodes2, 0);

	vector<double> bD = NeumannBoundaryCondition(meshFile, meshResult, tag, t, w);

	vector<double> x(nNodes2, 0);
	double xx = 0;

	double cD = 1e8;
	// Vyn�soben� vektoru koeficientem
	for (double& element : bD) {
		element *= cD;
	}
	double Q1 = 0.5 - sqrt(3.0 / 5.0) / 2.0;
	double Q2 = 0.5;
	double Q3 = 0.5 + sqrt(3.0 / 5.0) / 2.0;


	vector<vector<vector<double>>> basisFE_ref = {
	{{fiD(Q1), 0}, {fiE(Q1),0}, {fiF(Q1),0}, {0,fiD(Q1)} , {0,fiE(Q1)}, {0,fiF(Q1) }},
	{{fiD(Q2), 0}, {fiE(Q2),0}, {fiF(Q2),0}, {0,fiD(Q2)} , {0,fiE(Q2)}, {0,fiF(Q2) }},
	{{fiD(Q3), 0}, {fiE(Q3),0}, {fiF(Q3),0}, {0,fiD(Q3)} , {0,fiE(Q3)}, {0,fiF(Q3) }}
	};

	vector<double> w_g = { 5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0 };
	double nQuadrature = 3;
	int index_IJ = 0;

	for (int k = 0; k < vectorOfTags.size(); k++) {

		if (vectorOfTags[k] == tag) {

			int D = result1.matrix2[k][0] - 1;					// Ozna�en� vrchol� hran
			int E = result1.matrix2[k][1] - 1;
			int F = result1.matrix2[k][2] - 1;
			vector<int> edgeVertex = { D, E ,F };

			vector<double> cooD = result1.matrix3[D];		// Sou�adnice vrchol� hran
			vector<double> cooE = result1.matrix3[E];

			double edgeLenght = sqrt((cooD[0] - cooE[0]) * (cooD[0] - cooE[0]) + (cooD[1] - cooE[1]) * (cooD[1] - cooE[1]));
			vector<int> index = { D, E,F, D + nNodes, E + nNodes, F + nNodes };

			for (int i = 0; i < 6; i++) {
				for (int j = 0; j < 6; j++) {
					double k_val = 0;
					for (int l = 0; l < nQuadrature; l++) {
						vector<double> v = basisFE_ref[l][j];
						vector<double> u = basisFE_ref[l][i];

						//k_val = k_val + edgeLenght * W3dotprod(u, v);
						k_val = k_val + edgeLenght * w_g[l] * dotprod(u, v);

					};
					k_val = cD * k_val;
					I[index_IJ] = index[i];
					J[index_IJ] = index[j];
					VAL[index_IJ] = k_val;
					index_IJ++;
				};
			};
		};

	};
	SparseMatrixVectorResult assembledMatrix = AssembledMatrix(I, J, VAL, bD);



	return assembledMatrix;
}