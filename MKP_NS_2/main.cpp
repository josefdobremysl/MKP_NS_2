//12,12
#include "mesh_processor.h"
#include "Dirichlet.h"
#include "Neumann.h"
#include "Stiffnessmatrix.h"
#include "Massmatrix.h"
//#include "mesh_t.h"
#include "Aerodynamicforce.h"

#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <Eigen/IterativeLinearSolvers>

#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>	
#include "vtk.h"
#include "vtk_p.h"
#include <iomanip>

using namespace std;
using namespace Eigen;

int main() {

	// konstanty pro RK4
	//double  y = 0, dy = 0, ddy = 0, gama = 0.5, beta = 0.25, yn = 0, dyn = 0, ddyn = 0;

	int t_step = 0;
	int num_t_steps = 1000;
	//vector<double> vec_aero_F(num_t_steps, 0);
	double t = 0;
	double dt = 1./10.;

	string meshFile = "Obdelnik_O_1847_800.txt";

	MatrixVectorResult result = GetTriangles(meshFile);
	int nTri = result.integerResult;		// poèet trojúhelníkù
	int nNodes = result.matrix3.size();			// poèet bodù

	vector<int> MMM(nNodes, 0);
	for (int k = 0; k < nTri; k++) {
		int A = result.matrix1[k][0] - 1;					// Oznaèení vrcholù trojúhelníku k
		int B = result.matrix1[k][1] - 1;
		int C = result.matrix1[k][2] - 1;
		vector<int> verTri = { A, B, C };
		for (int i = 0; i < 3; i++) {
			MMM[verTri[i]] = 1;
		}
	}

	int nVert = 0;

	for (int k = 0; k < nNodes; k++) {
		nVert += MMM[k];
	}
	int xsize = 2 * nNodes + nVert;

	VectorXd x = VectorXd::Zero(xsize);
	VectorXd x_k = VectorXd::Zero(xsize);
	//x.segment(0, nNodes).setOnes();
	//x_k.segment(0, nNodes).setOnes();

	//VectorXd displac_p = VectorXd::Zero(xsize);
	//VectorXd displac_n = VectorXd::Zero(xsize);
	std::vector<double> wvec(xsize, 0.0);
	//Eigen::VectorXd wvec = Eigen::VectorXd::Zero(xsize);
	//Eigen wvec = Eigen::Zero(xsize);


	int nNodes2 = 2 * nNodes;
	vector<vector<double>>NodesM = result.matrix3;
	VectorXd xPrev; // Uchování předchozího řešení
	bool converged = false;
	int iteration = 0;




	for (int j = 0; j < num_t_steps; j++) {
		if (j > 50) {
			dt=0.05;
		}
		//for (int i = 0; i < 1; i++) {
		//MatrixVectorResult result1 = GetTriangles(meshFile);
		//DisplacementResult displ_res = computeDisplacement(meshFile, result1, t, yn);
		//MatrixVectorResult meshResult = MeshT(meshFile, result, displ_res);
		//displac_n.head(displ_res.displacements.size()) = displ_res.displacements; // Zkopírování hodnot z displ_res.displacements

		//w = (displac_n - displac_p) / dt;
		//vector<double> wvec(w.data(), w.data() + w.size());
		//displac_p = displac_n;


		SparseMatrixVectorResult mass_matrix = MassMatrix(meshFile, result);
		SparseMatrix<double> M = mass_matrix.sparsematrix;



		// Apply Dirichlet boundary conditions
		int tags[] = { 100,400, 300, 500, 600 };
		Eigen::SparseMatrix<double> sparsematrixDirichlet(xsize, xsize);
		vector<double> combinedBoundaryVector(xsize, 0.0);

		for (int tag : tags) {
			SparseMatrixVectorResult Dirichlet = DirichletBoundaryCondition(meshFile, result, tag, t, wvec);
			sparsematrixDirichlet += Dirichlet.sparsematrix;
			const auto& bD = Dirichlet.right_side_vector;
			for (size_t i = 0; i < bD.size(); ++i) {
				combinedBoundaryVector[i] += bD[i];
			}
		}


		while (!converged) {

			iteration++;

			SparseMatrixVectorResult stiffness_matrix = StiffnessMatrix(meshFile, result, x_k, t);
			SparseMatrix<double> K = stiffness_matrix.sparsematrix;
			vector<double> right_side_vector0 = stiffness_matrix.right_side_vector;

			/////////////////////////////////////////////////// Crank - Nikolson
			//vector<double> right_side_vector(xsize, 0);
			//Eigen::VectorXd dtMv = ((M / dt - K / 2) * x);
			//for (int i = 0; i < xsize; i++) {
			//	//right_side_vector[i] = right_side_vector0[i] + combinedBoundaryVector[i] + dtMv[i];
			//	right_side_vector[i] = right_side_vector0[i] + combinedBoundaryVector[i];

			//}
			//SparseMatrix<double> A = M / dt + K / 2 + sparsematrixDirichlet;


			///////////////////////////////////////////////// implicitni Euler
			vector<double> right_side_vector(xsize, 0);
			VectorXd dtMv = ((M / dt) * x);
			for (int i = 0; i < xsize; i++) {
				right_side_vector[i] = right_side_vector0[i] + combinedBoundaryVector[i] +dtMv[i];
				//right_side_vector[i] = right_side_vector0[i] + combinedBoundaryVector[i];
			}
			SparseMatrix<double> A = M / dt + sparsematrixDirichlet + K;
			//SparseMatrix<double> A = sparsematrixDirichlet + K;
			/////////////////////////////////////////////////////////////

					// Pravá strana soustavy b
			VectorXd b = Map<VectorXd>(right_side_vector.data(), right_side_vector.size());

			//////////////////////  ŘEŠENÍ SOUSTAVY Konjugovane gradienty - Eigen\\\\\\\\\\\\\\\\\\\\\\\\\\
				//// Inicializace øešení
			SparseLU<SparseMatrix<double>> solver;
			solver.compute(A);
			xPrev = x_k; // Uložení předchozího řešení
			x_k = solver.solve(b);

			//std::cout << MatrixXd(K) << std::endl;

			// Tisk aktuální normy rozdílu
			double diffNorm = (x_k - xPrev).norm();
			std::cout << "Iterace " << iteration << ": Norma rozdílu = " << diffNorm << std::endl;


			// Kontrola konvergence
			converged = (diffNorm < 0.005);
			if (converged == true) {
				x = x_k;
			}


		}
		converged = false;

		size_t dotPos1 = meshFile.find_last_of(".");
		string File_Aero_F = "Aero_Force_"+ meshFile.substr(0, dotPos1) + ".txt";
		vector<double> xvec(x.data(), x.data() + x.size());
		vector<double> FF5 = AerodynamicForce(meshFile, result, 500, t, xvec);
		vector<double> FF6 = AerodynamicForce(meshFile, result, 600, t, xvec);
		double FF = FF5[1] + FF6[1];
		std::ofstream file(File_Aero_F, std::ios::app);
		
		if (!file) {
			std::cerr << "Nelze otevrit soubor pro zapis.\n";
		}
		file << j << " " << FF5[1] << " " << FF6[1] << " " << FF << "\n";
		
		file.close();
		std::cout << "Data byla zapsana do souboru \n";

		//double pi = 3.141592;
		//double Ms = 2;
		//double Ur = 3;
		//double CL = 30;
		//double b = CL / 2 / Ms;
		//double a = 4 * pi * pi / Ur / Ur;

		//vec_aero_F[t_step] = FF;
		/////////////   Newmarkova metoda   ///////////////

		//ddyn = b * FF - a * y;
		//yn = 0;//y + dt * dy + dt * dt / 2 * ((1 - 2 * beta) * ddy + 2 * beta * ddyn);
		//dyn = dy + (1 - gama) * dt * ddy + gama * dt * ddyn;
		//dy = dyn;
		//y = yn;
		//y = 0;
		//}

			///////////////////////////////		VYPIS VTK		\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

			//MatrixVectorResult result1 = GetTriangles(meshFile);
			//DisplacementResult displ_res = computeDisplacement(meshFile, result1, t, yn);
			//MatrixVectorResult meshResult = MeshT(meshFile, result, displ_res);
		string meshFile2 = meshFile;

		// Odvození jmena vystupniho souboru
		size_t dotPos = meshFile2.find_last_of(".");

		if (dotPos != std::string::npos && dotPos != meshFile2.size() - 1) {
			// Nahrazení části řetězce začínající po tečce novým řetězcem
			meshFile2 = meshFile.substr(0, dotPos) + "_" + std::to_string(j) + ".vtk";
		}
		else {
			meshFile2 += ".vtk";
		}
		VTKfile(meshFile2, meshFile, result, x);

		vector<double> p(x.size() - nNodes2, 0);
		int ii = 0;
		for (int i = nNodes2; i < x.size(); i++) {
			p[ii] = x[i];
			ii++;
		}
		string meshFilePress = meshFile.substr(0, dotPos) + "_" + std::to_string(j) + "pressure.vtk";


		VectorXd pp = Map<VectorXd>(p.data(), p.size());
		VTKfile_p(meshFilePress, meshFile, result, pp);

		t += dt;
		t_step++;
		//////////////////////////////////////////////////////////////////////////////////////////////


	}
}