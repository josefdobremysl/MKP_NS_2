#ifndef FUNCTION_H
#define FUNCTION_H

#include "mesh_processor.h"

#include <string>
#include <vector>
#include <tuple>
#include <Eigen/Sparse>
#include <armadillo>

using namespace Eigen;
using namespace std;
using namespace arma;

double Function(double& X, double& Y);
double FunctionN1(double& X, double& Y, const int& tag, vector<double>& n, vector<double> w);
double FunctionN2(double& X, double& Y, const int& tag, vector<double>& n, vector<double> w);
double fi1(double& x);
double fi2(double& x);
double fiD(double& x);
double fiE(double& x);
double fiF(double& x);


double dotprod(vector<double>& u, vector<double>& v);
double W3dotprod(vector<double>& u, vector<double>& v);

vector<double> FunctionV(double& X, double& Y);

mat gradFEx_ref_Assem(vector<int>& ind);
mat gradFEy_ref_Assem(vector<int>& ind);
mat gradFEx_ref_Assem3(vector<int>& ind);
mat gradFEy_ref_Assem3(vector<int>& ind);
mat basis_Assem(vector<int>& idx);


double gradFEx_ref2(const int& i, const int& j);
double gradFEy_ref2(const int& i, const int& j);
double basis_2(const int& i, const int& j);
double gradFEx_ref3(const int& i, const int& j);
double gradFEy_ref3(const int& i, const int& j);

#endif // FUNCTION