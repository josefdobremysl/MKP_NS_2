#include "Massmatrix.h"

#include "mesh_processor.h"

#include "function.h"
#include "Assemblematrix.h"

#include <Eigen/Sparse>
#include <armadillo>
#include <omp.h>


#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <tuple>
#include <numeric>
#include <chrono>	
using namespace Eigen;
using namespace std;
using namespace arma;



SparseMatrixVectorResult MassMatrix(const string& meshFile, MatrixVectorResult& meshResult) {

	MatrixVectorResult result = meshResult;

	int nTri = result.integerResult;		// po�et troj�heln�k�
	int nNodes = result.matrix3.size();			// po�et bod�

	vector<int> M(nNodes, 0);
	for (int k = 0; k < nTri; k++) {
		int A = result.matrix1[k][0] - 1;					// Ozna�en� vrchol� troj�heln�ku k
		int B = result.matrix1[k][1] - 1;
		int C = result.matrix1[k][2] - 1;
		vector<int> verTri = { A, B, C };
		for (int i = 0; i < 3; i++) {
			M[verTri[i]] = 1;
		}
	}

	int nVert = 0;

	for (int k = 0; k < nNodes; k++) {
		nVert += M[k];
	}
	//cout << "index_pressure";
	int trans_i = 0;
	vector<int> index_pressure(nNodes, 0);
	for (int i = 0; i < nNodes; i++) {
		if (M[i] != 0) {
			index_pressure[i] = trans_i;
			trans_i++;
		}

	}
	cout << endl;

	//cout << "nNodes  " << nNodes << endl;
	int nNodes2 = 2 * nNodes;
	int nNodesIJ = nNodes2 * 12 * 12 * 4;

	vector<int> I(nNodesIJ, 0);
	vector<int> J(nNodesIJ, 0);
	vector<double> VAL(nNodesIJ, 0);
	vector<double> right_side_vector(nNodes2 + nVert, 0);		// vektor prav� strany


	// Vytvo�en� matice basisFE_ref
	mat basisFE_ref = { // bazove f-ce v x-smeru a pak baz. fce v y-smeru
		{0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
		{0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
		{0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
		{0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
		{0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
		{0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
		{0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
		{0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
		{0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
		{0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0},
		{0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0},
		{0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0}
	};


	double nQuadrature = 3;
	double w = 1 / nQuadrature;
	int index_IJ = 0;


	for (int k = 0; k < nTri; k++) {
		int A = result.matrix1[k][0] - 1;					// Ozna�en� vrchol� troj�heln�ku k
		int B = result.matrix1[k][1] - 1;
		int C = result.matrix1[k][2] - 1;
		int D = result.matrix1[k][3] - 1;					// Ozna�en� stredu stran troj�heln�ku k
		int E = result.matrix1[k][4] - 1;
		int F = result.matrix1[k][5] - 1;


		vector<int> verTri = { A, B, C, D, E, F };

		vector<double> cooA = result.matrix3[A];		// Sou�adnice vrchol� troj�heln�ku k
		vector<double> cooB = result.matrix3[B];
		vector<double> cooC = result.matrix3[C];

		mat CoQ = {
			{ cooA[0],cooB[0],cooC[0],(cooA[0] + cooB[0]) / 2, (cooB[0] + cooC[0]) / 2, (cooC[0] + cooA[0]) / 2},
			{ cooA[1],cooB[1],cooC[1],(cooA[1] + cooB[1]) / 2, (cooB[1] + cooC[1]) / 2, (cooC[1] + cooA[1]) / 2},
		};

		mat matA = {
			{ cooB[0] - cooA[0], cooC[0] - cooA[0] },
			{ cooB[1] - cooA[1], cooC[1] - cooA[1] },
		};

		double detA = matA(0, 0) * matA(1, 1) - matA(0, 1) * matA(1, 0);

		mat invA = {
			{  matA(1,1) / detA, -matA(0,1) / detA },
			{ -matA(1,0) / detA,  matA(0,0) / detA  },
		};


		vector<int> index = { A, B, C, D, E, F, nNodes + A, nNodes + B, nNodes + C, nNodes + D, nNodes + E, nNodes + F };

		for (int i = 0; i < 12; ++i) {

			for (int j = 0; j < 12; ++j) {
				double k_val = 0;
				for (int l = 0; l < 12; ++l) {

					k_val = k_val + detA / 2 * w * basisFE_ref(l, i) * basisFE_ref(l, j);

				};

				// zapis vysledneho soucinu vektoru-tenzoru na spravne misto do tripletu 
				I[index_IJ] = index[i];
				J[index_IJ] = index[j];
				VAL[index_IJ] = k_val;
				index_IJ++;
			}

		};
	};

	SparseMatrixVectorResult assembledMatrix = AssembledMatrix(I, J, VAL, right_side_vector);
	SparseMatrix<double> K = assembledMatrix.sparsematrix;

	return assembledMatrix;
}