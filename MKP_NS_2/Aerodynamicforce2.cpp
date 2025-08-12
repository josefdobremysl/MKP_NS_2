#include "mesh_processor.h"
#include "Aerodynamicforce2.h"

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



SparseMatrixVectorResult AerodynamicForce2(const string& meshFile, MatrixVectorResult& meshResult, Eigen::VectorXd x, const double& t) {

	double Re = 800.;
	double mu = 1. / Re;


	MatrixVectorResult result = meshResult;

	int nTri = result.integerResult;		// po�et troj�heln�k�
	int nNodes = result.matrix3.size();			// pocet bodu

	vector<int> M(nNodes, 0);
	for (int k = 0; k < nTri; k++) {
		int A = result.matrix1[k][0] - 1;					// Oznaceni vrcholu trojuhelniku k
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

	double F_d = 0.0;
	double F_n = 0.0;


	mat basis_q_ref = {
	{0.5,  0.5,  0.0},
	{0.0,  0.5,  0.5},
	{0.5,  0.0,  0.5},
	};



	//vector<int> idx1 = {0,1,2, 3, 4, 5,6,7,8, 9, 10, 11,12,13 };
	vector<int> idx1 = { 0,1,2, 3, 4, 5,6,7 };
	//vector<int> idx1 = { 0,1,2, 3, 4, 5};
	// Vytvo�en� matice basisFE_ref
	mat basisFE_ref = basis_Assem(idx1); // bazove f-ce v x-smeru a pak baz. fce v y-smeru

	// PD podle X baz. fci v kvadr. uzlech
	mat gradFEx_ref = gradFEx_ref_Assem(idx1);
	// PD podle Y baz. fci v kvadr. uzlech
	mat gradFEy_ref = gradFEy_ref_Assem(idx1);

	vector<int> idx3 = { 0, 1, 2,3,4,5 };
	// PD podle X baz. fci v kvadr. uzlech
	mat gradFEx_ref_3 = gradFEx_ref_Assem3(idx3);
	// PD podle Y baz. fci v kvadr. uzlech
	mat gradFEy_ref_3 = gradFEy_ref_Assem3(idx3);

	//vector<int> idx2 = { 0, 1, 2, 6, 7, 8 };
	//// PD podle X baz. fci v kvadr. uzlech
	//mat gradFEx_refV = gradFEx_ref_Assem(idx2);
	//// PD podle Y baz. fci v kvadr. uzlech
	//mat gradFEy_refV = gradFEx_ref_Assem(idx2);
	double nQuadrature = 3;
	double nQuadrature2 = 4;

	double w3 = 1. / 3.;
	vector<double> w7 = { 0.05,0.05,0.05,8. / 60.,8. / 60.,8. / 60.,27. / 60. };
	vector<double> w4 = { 25. / 96.,25. / 96.,25. / 96.,-9. / 32. };

	int index_IJ = 0;


	for (int k = 0; k < nTri; k++) {
		int A = result.matrix1[k][0] - 1;					// Ozna�en� vrchol� troj�heln�ku k
		int B = result.matrix1[k][1] - 1;
		int C = result.matrix1[k][2] - 1;
		int D = result.matrix1[k][3] - 1;					// Ozna�en� stredu stran troj�heln�ku k
		int E = result.matrix1[k][4] - 1;
		int F = result.matrix1[k][5] - 1;

		double xG_x = -1. / 9. * x(A) - 1. / 9. * x(B) - 1. / 9. * x(C) + 4. / 9. * x(D) + 4. / 9. * x(E) + 4. / 9. * x(F);
		double xG_y = -1. / 9. * x(A + nNodes) - 1. / 9. * x(B + nNodes) - 1. / 9. * x(C + nNodes) + 4. / 9. * x(D + nNodes) + 4. / 9. * x(E + nNodes) + 4. / 9. * x(F + nNodes);

		//////////////
		double vQ1x = 2. / 9. * x(A) - 1. / 9. * x(B) - 1. / 9. * x(C) + 4. / 9. * x(D) + 1. / 9. * x(E) + 4. / 9. * x(F);
		double vQ2x = -1. / 9. * x(A) + 2. / 9. * x(B) - 1. / 9. * x(C) + 4. / 9. * x(D) + 4. / 9. * x(E) + 1. / 9. * x(F);
		double vQ3x = -1. / 9. * x(A) - 1. / 9. * x(B) + 2. / 9. * x(C) + 1. / 9. * x(D) + 4. / 9. * x(E) + 4. / 9. * x(F);
		double vQ4x = -1. / 9. * (x(A) + x(B) + x(C)) + 4. / 9. * (x(D) + x(E) + x(F));

		double vQ1y = 2. / 9. * x(A + nNodes) - 1. / 9. * x(B + nNodes) - 1. / 9. * x(C + nNodes) + 4. / 9. * x(D + nNodes) + 1. / 9. * x(E + nNodes) + 4. / 9. * x(F + nNodes);
		double vQ2y = -1. / 9. * x(A + nNodes) + 2. / 9. * x(B + nNodes) - 1. / 9. * x(C + nNodes) + 4. / 9. * x(D + nNodes) + 4. / 9. * x(E + nNodes) + 1. / 9. * x(F + nNodes);
		double vQ3y = -1. / 9. * x(A + nNodes) - 1. / 9. * x(B + nNodes) + 2. / 9. * x(C + nNodes) + 1. / 9. * x(D + nNodes) + 4. / 9. * x(E + nNodes) + 4. / 9. * x(F + nNodes);
		double vQ4y = -1. / 9. * (x(A + nNodes) + x(B + nNodes) + x(C + nNodes)) + 4. / 9. * (x(D + nNodes) + x(E + nNodes) + x(F + nNodes));

		mat vQ = { vQ1x, vQ2x, vQ3x, vQ4x, vQ1y, vQ2y, vQ3y, vQ4y };



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

		mat gradFEx2 = invA(0, 0) * gradFEx_ref + invA(1, 0) * gradFEy_ref;
		mat gradFEy2 = invA(0, 1) * gradFEx_ref + invA(1, 1) * gradFEy_ref;

		mat gradFEx3 = invA(0, 0) * gradFEx_ref_3 + invA(1, 0) * gradFEy_ref_3;
		mat gradFEy3 = invA(0, 1) * gradFEx_ref_3 + invA(1, 1) * gradFEy_ref_3;


		vector<int> index = { A, B, C, D, E, F, nNodes + A, nNodes + B, nNodes + C, nNodes + D, nNodes + E, nNodes + F };

		for (int i = 0; i < 12; ++i) {

			for (int j = 0; j < 12; ++j) {
				double k_val = 0;
				double D = 0;
				if ((i < 6 || j < 6) ^ (i > 5 || j > 5)) {
					D = 1;
				}
				for (int l = 0; l < nQuadrature2; ++l) {

					vector<double> gradFI_i = { gradFEx2(l,i), gradFEy2(l,i), gradFEx2(l + nQuadrature2,i), gradFEy2(nQuadrature2 + l,i) };
					vector<double> gradFI_j = { gradFEx2(l,j), gradFEy2(l,j), gradFEx2(l + nQuadrature2,j), gradFEy2(nQuadrature2 + l,j) };

					vector<double> vQl = { vQ(l),vQ(l + nQuadrature2),vQ(l),vQ(l + nQuadrature2) };

					double konvekce_i = dotprod(vQl, gradFI_j);


					k_val = k_val + detA / 2. * w4[l] * (mu * dotprod(gradFI_i, gradFI_j) + konvekce_i * basis_2(l, i) * D - pQ);

				};
				F_d = F_d + k_val;				
				// zapis vysledneho soucinu vektoru-tenzoru na spravne misto do tripletu 

				index_IJ++;
			}

		};
	};



	return assembledMatrix;
}