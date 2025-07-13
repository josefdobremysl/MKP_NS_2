#include "mesh_processor.h"
#include "Dirichlet.h"
#include "Neumann.h"
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

vector<double> NeumannBoundaryCondition(const string& meshFile, MatrixVectorResult& meshResult, const int& tag, const double& t, vector<double> w) {

	MatrixVectorResult loadedMesh = meshResult;
	int nTri = loadedMesh.integerResult;		// po�et troj�heln�k�
	int nNod = loadedMesh.matrix3.size();
	vector<int> M(nNod, 0);
	for (int k = 0; k < nTri; k++) {
		int A = loadedMesh.matrix1[k][0] - 1;					// Ozna�en� vrchol� troj�heln�ku k
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
	vector<double> bNinPoints(nNodes, 0);


	for (int i = 0; i < size; i++) {
		if (vectorOfTags[i] == tag) {
			int D = lineVertex[i][0] - 1;
			int E = lineVertex[i][1] - 1;
			int F = lineVertex[i][2] - 1;

			vector<double> cooD = loadedMesh.matrix3[D];		// Sou�adnice vrchol� troj�heln�ku k
			vector<double> cooE = loadedMesh.matrix3[E];

			vector<double> wD = { w[D],w[D + N] };
			vector<double> wE = { w[E],w[E + N] };

			vector<double> n = { cooE[1] - cooD[1], cooD[0] - cooE[0] };


			//////////////////////////////////////////////////////////////////////////////

			double Q1 = 0.5 - sqrt(3.0 / 5.0) / 2.0;
			double Q2 = 0.5;
			double Q3 = 0.5 + sqrt(3.0 / 5.0) / 2.0;


			vector<double> cooF1 = { cooD[0] + (cooE[0] - cooD[0]) * Q1, cooD[1] + (cooE[1] - cooD[1]) * Q1 };
			vector<double> cooF2 = { cooD[0] + (cooE[0] - cooD[0]) * Q2, cooD[1] + (cooE[1] - cooD[1]) * Q2 };
			vector<double> cooF3 = { cooD[0] + (cooE[0] - cooD[0]) * Q3, cooD[1] + (cooE[1] - cooD[1]) * Q3 };

			vector<double> wF1 = { wD[0] + (wE[0] - wD[0]) * Q1, wD[1] + (wE[1] - wD[1]) * Q1 };
			vector<double> wF2 = { wD[0] + (wE[0] - wD[0]) * Q2, wD[1] + (wE[1] - wD[1]) * Q2 };
			vector<double> wF3 = { wD[0] + (wE[0] - wD[0]) * Q3, wD[1] + (wE[1] - wD[1]) * Q3 };

			double edgeLength = sqrt(pow(cooD[0] - cooE[0], 2) + pow(cooD[1] - cooE[1], 2));

			vector<double> fN1 = { FunctionN1(cooF1[0], cooF1[1], tag, n, wF1), FunctionN1(cooF2[0], cooF2[1], tag, n, wF2) ,FunctionN1(cooF3[0], cooF3[1], tag, n, wF3) };
			vector<double> fN2 = { FunctionN2(cooF1[0], cooF1[1], tag, n, wF1), FunctionN2(cooF2[0], cooF2[1], tag, n, wF2) ,FunctionN2(cooF3[0], cooF3[1], tag, n, wF3) };

			vector<double> fi_D_F1 = { fiD(Q1),fiD(Q2),fiD(Q3) };
			vector<double> fi_D_F2 = { fiD(Q1),fiD(Q2),fiD(Q3) };
			vector<double> fi_E_F1 = { fiE(Q1),fiE(Q2),fiE(Q3) };
			vector<double> fi_E_F2 = { fiE(Q1),fiE(Q2),fiE(Q3) };
			vector<double> fi_F_F1 = { fiF(Q1),fiF(Q2),fiF(Q3) };
			vector<double> fi_F_F2 = { fiF(Q1),fiF(Q2),fiF(Q3) };

			double w = 1.;

			bNinPoints[D] += edgeLength * w * W3dotprod(fN1, fi_D_F1);
			bNinPoints[E] += edgeLength * w * W3dotprod(fN1, fi_E_F1);
			bNinPoints[F] += edgeLength * w * W3dotprod(fN1, fi_F_F1);
			bNinPoints[D + N] += edgeLength * w * W3dotprod(fN2, fi_D_F2);
			bNinPoints[E + N] += edgeLength * w * W3dotprod(fN2, fi_E_F2);
			bNinPoints[F + N] += edgeLength * w * W3dotprod(fN2, fi_F_F2);


		}
	}


	return bNinPoints;

};