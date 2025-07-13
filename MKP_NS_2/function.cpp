#include "mesh_processor.h"
#include "function.h"
//#include "print.h"
//#include "quicksort.h"

#include <Eigen/Sparse>
#include <armadillo>

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <tuple>


//using namespace Eigen;
using namespace std;
using namespace arma;

const double pi = 3.14159265358979323846264338327950288;

double f1;
double f2;

////////////////////////////////////////////////////////////////////////////////////
///////////    Case1
vector<double> FunctionV(double& X, double& Y) {
	vector<double> fV(2, 0);
	fV[0] = 0;
	fV[1] = 0;
	return fV;

}



double FunctionN1(double& X, double& Y, const int& tag, vector<double>& n, vector<double> w) {

	double n_norm = (sqrt(n[0] * n[0] + n[1] * n[1]));
	vec n_jednotkovy = { n[0] / n_norm,n[1] / n_norm };

	if (tag == 100) {
		f1 = 0;
		//f1 = 1 - ((Y-1) / 0.5 - 1) * ((Y - 1) / 0.5 - 1);

	}
	if (tag == 200) {
		//f1 = 1 - (Y / 0.5 - 1) * (Y / 0.5 - 1);
		f1 = 0;// 1 - (Y / 0.5 - 1) * (Y / 0.5 - 1);// 1 * sin(pi * Y);
	}
	if (tag == 300) {
		f1 = 0;
	}
	if (tag == 400) {
		f1 = 1 - (Y / 0.5 - 1) * (Y / 0.5 - 1);
	}
	if (tag == 500) {
		f1 = 0;// w[0];
		//cout << "w0 " << w[0] << endl;
	}
	if (tag == 600) {
		f1 = 0;//w[0];// -1 * sin(pi * Y);
		//cout << "w0 " << w[0] << endl;
	}

	return f1;
}

double FunctionN2(double& X, double& Y, const int& tag, vector<double>& n, vector<double> w) {

	double n_norm = (sqrt(n[0] * n[0] + n[1] * n[1]));
	vec n_jednotkovy = { n[0] / n_norm,n[1] / n_norm };

	if (tag == 100) {
		f2 = 0;
	}
	if (tag == 200) {
		f2 = 0;
	}
	if (tag == 300) {
		f2 = 0;
	}
	if (tag == 400) {
		f2 = 0;
	}
	if (tag == 500) {
		f2 = 0;// w[1];
		//cout << "w1 " << w[1]*1000 << endl;
	}
	if (tag == 600) {
		f2 = 0; // w[1];
		//cout << "w1 " << w[1] * 1000 << endl;
	}
	return f2;
}



//bases functions
double fi1(double& x) {
	double fi1 = 1 - x;
	return fi1;
}
double fi2(double& x) {
	double fi2 = x;
	return fi2;
}
double fiD(double& x) {
	double fiD = 1 - 3 * x + 2 * x * x;
	return fiD;
}
double fiE(double& x) {
	double fiE = 2 * x * x - x;
	return fiE;
}
double fiF(double& x) {
	double fiF = 4 * x - 4 * x * x;
	return fiF;
}



double dotprod(vector<double>& u, vector<double>& v) {
	double result = 0;
	for (int i = 0; i < u.size(); ++i) {
		result += u[i] * v[i];
	}
	return result;
}

double W3dotprod(vector<double>& u, vector<double>& v) {
	double result1 = 0;
	result1 = 5.0 / 18.0 * u[0] * v[0] + 8.0 / 18.0 * u[1] * v[1] + 5.0 / 18.0 * u[2] * v[2];
	return result1;
}

//double basis_2(const int& i, const int& j) {
//	double x = 0;
//	double y = 0;
//	double bas = 0;
//	if (i == 0 || i == 7) { x = 0.;   y = 0.; }
//	if (i == 1 || i == 8) { x = 1.;   y = 0.; }
//	if (i == 2 || i == 9) { x = 0.;    y = 1.; }
//	if (i == 3 || i == 10) { x = 0.5; y = 0.; }
//	if (i == 4 || i == 11) { x = 0.5;   y = 0.5; }
//	if (i == 5 || i == 12) { x = 0.;   y = 0.5; }
//	if (i == 6 || i == 13) { x = 0.3;   y = 0.3; }
//
//	if (j == 0 && i < 7 || j == 6 && i > 6) { bas = (1 - x - y) * (1 - 2 * x - 2 * y); }
//	if (j == 1 && i < 7 || j == 7 && i > 6) { bas = x * (2 * x - 1); }
//	if (j == 2 && i < 7 || j == 8 && i > 6) { bas = y * (2 * y - 1); }
//	if (j == 3 && i < 7 || j == 9 && i > 6) { bas = 4 * x * (1 - x - y); }
//	if (j == 4 && i < 7 || j == 10 && i > 6) { bas = 4 * x * y; }
//	if (j == 5 && i < 7 || j == 11 && i > 6) { bas = 4 * y * (1 - x - y); }
//	return bas;
//}
//
//mat basis_Assem(vector<int>& idx) {
//	mat A0(14, 12);
//
//	for (int i : idx) {
//		for (int j = 0; j < 12; ++j) {
//			A0(i, j) = basis_2(i, j);
//		}
//	}
//	return A0;
//}
//
//mat gradFEx_ref_Assem(vector<int>& idx) {
//	mat A1(14, 12);
//
//	for (int i : idx) {
//		for (int j = 0; j < 12; ++j) {
//			A1(i, j) = gradFEx_ref2(i, j);
//		}
//	}
//	return A1;
//}
//mat gradFEy_ref_Assem(vector<int>& idx) {
//	mat A2(14, 12);
//
//	for (int i : idx) {
//		for (int j = 0; j < 12; ++j) {
//			A2(i, j) = gradFEy_ref2(i, j);
//		}
//	}
//	return A2;
//}
//
mat gradFEx_ref_Assem3(vector<int>& idx) {
	mat A1(6, 12);

	for (int i : idx) {
		for (int j = 0; j < 12; ++j) {
			A1(i, j) = gradFEx_ref3(i, j);
		}
	}
	return A1;
}
mat gradFEy_ref_Assem3(vector<int>& idx) {
	mat A2(6, 12);

	for (int i : idx) {
		for (int j = 0; j < 12; ++j) {
			A2(i, j) = gradFEy_ref3(i, j);
		}
	}
	return A2;
}
//
//double gradFEx_ref2(const int& i, const int& j) {
//	double x = 0;
//	double y = 0;
//	double grx = 0;
//	if (i == 0 || i == 7) { x = 0.;   y = 0.; }
//	if (i == 1 || i == 8) { x = 1.;   y = 0.; }
//	if (i == 2 || i == 9) { x = 0.;    y = 1.; }
//	if (i == 3 || i == 10) { x = 0.5; y = 0.; }
//	if (i == 4 || i == 11) { x = 0.5;   y = 0.5; }
//	if (i == 5 || i == 12) { x = 0.;   y = 0.5; }
//	if (i == 6 || i == 13) { x = 0.3;   y = 0.3; }
//
//	if (j == 0 && i < 7 || j == 6 && i > 6) { grx = 4 * x + 4 * y - 3; }
//	if (j == 1 && i < 7 || j == 7 && i > 6) { grx = 4 * x - 1; }
//	if (j == 2 && i < 7 || j == 8 && i > 6) { grx = 0; }
//	if (j == 3 && i < 7 || j == 9 && i > 6) { grx = 4 - 8 * x - 4 * y; }
//	if (j == 4 && i < 7 || j == 10 && i > 6) { grx = 4 * y; }
//	if (j == 5 && i < 7 || j == 11 && i > 6) { grx = -4 * y; }
//	return grx;
//}
//
//double gradFEy_ref2(const int& i, const int& j) {
//	double x = 0;
//	double y = 0;
//	double gry = 0;
//	if (i == 0 || i == 7) { x = 0.;   y = 0.; }
//	if (i == 1 || i == 8) { x = 1.;   y = 0.; }
//	if (i == 2 || i == 9) { x = 0.;    y = 1.; }
//	if (i == 3 || i == 10) { x = 0.5; y = 0.; }
//	if (i == 4 || i == 11) { x = 0.5;   y = 0.5; }
//	if (i == 5 || i == 12) { x = 0.;   y = 0.5; }
//	if (i == 6 || i == 13) { x = 0.3;   y = 0.3; }
//
//	if (j == 0 && i < 7 || j == 6 && i > 6) { gry = 4 * x + 4 * y - 3; }
//	if (j == 1 && i < 7 || j == 7 && i > 6) { gry = 0; }
//	if (j == 2 && i < 7 || j == 8 && i > 6) { gry = 4 * y - 1; }
//	if (j == 3 && i < 7 || j == 9 && i > 6) { gry = -4 * x; }
//	if (j == 4 && i < 7 || j == 10 && i > 6) { gry = 4 * x; }
//	if (j == 5 && i < 7 || j == 11 && i > 6) { gry = 4 - 4 * x - 8 * y; }
//	return gry;
//}
double gradFEx_ref3(const int& i, const int& j) {
	double x = 0.;
	double y = 0.;
	double grx = 0.;
	if (i == 0 || i == 3) { x = 0.5;   y = 0.; }
	if (i == 1 || i == 4) { x = 0.5;   y = 0.5; }
	if (i == 2 || i == 5) { x = 0.0;   y = 0.5; }
	

	if (j == 0 && i < 3 || j == 6  && i > 2) { grx = 4 * x + 4 * y - 3; }
	if (j == 1 && i < 3 || j == 7  && i > 2) { grx = 4 * x - 1; }
	if (j == 2 && i < 3 || j == 8  && i > 2) { grx = 0; }
	if (j == 3 && i < 3 || j == 9  && i > 2) { grx = 4 - 8 * x - 4 * y; }
	if (j == 4 && i < 3 || j == 10 && i > 2) { grx = 4 * y; }
	if (j == 5 && i < 3 || j == 11 && i > 2) { grx = -4 * y; }
	return grx;
}

double gradFEy_ref3(const int& i, const int& j) {
	double x = 0.;
	double y = 0.;
	double gry = 0.;
	if (i == 0 || i == 3) { x = 0.5;   y = 0.; }
	if (i == 1 || i == 4) { x = 0.5;   y = 0.5; }
	if (i == 2 || i == 5) { x = 0.0;   y = 0.5; }


	if (j == 0 && i < 3 || j == 6 && i > 2) { gry = 4 * x + 4 * y - 3; }
	if (j == 1 && i < 3 || j == 7 && i > 2) { gry = 0; }
	if (j == 2 && i < 3 || j == 8 && i > 2) { gry = 4 * y - 1; }
	if (j == 3 && i < 3 || j == 9 && i > 2) { gry = -4 * x; }
	if (j == 4 && i < 3 || j == 10 && i > 2) { gry = 4 * x; }
	if (j == 5 && i < 3 || j == 11 && i > 2) { gry = 4 - 4 * x - 8 * y; }
	return gry;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
double basis_2(const int& i, const int& j) {
	double x = 0.;
	double y = 0.;
	double bas = 0.;
	if (i == 0 || i == 4) { x = 1./6.;   y = 1. / 6.; }
	if (i == 1 || i == 5) { x = 2./3.;   y = 1. / 6.; }
	if (i == 2 || i == 6) { x = 1./6.;   y = 2./3.;   }
	if (i == 3 || i == 7) { x = 1./3.;   y = 1. / 3.; }
	

	if (j == 0 && i < 4 || j == 6 && i > 3) { bas = (1 - x - y) * (1 - 2 * x - 2 * y); }
	if (j == 1 && i < 4 || j == 7 && i > 3) { bas = x * (2 * x - 1); }
	if (j == 2 && i < 4 || j == 8 && i > 3) { bas = y * (2 * y - 1); }
	if (j == 3 && i < 4 || j == 9 && i > 3) { bas = 4 * x * (1 - x - y); }
	if (j == 4 && i < 4 || j == 10 && i > 3) { bas = 4 * x * y; }
	if (j == 5 && i < 4 || j == 11 && i > 3) { bas = 4 * y * (1 - x - y); }
	return bas;
}

mat basis_Assem(vector<int>& idx) {
	mat A0(8, 12);

	for (int i : idx) {
		for (int j = 0; j < 12; ++j) {
			A0(i, j) = basis_2(i, j);
		}
	}
	return A0;
}

mat gradFEx_ref_Assem(vector<int>& idx) {
	mat A1(8, 12);

	for (int i : idx) {
		for (int j = 0; j < 12; ++j) {
			A1(i, j) = gradFEx_ref2(i, j);
		}
	}
	return A1;
}
mat gradFEy_ref_Assem(vector<int>& idx) {
	mat A2(8, 12);

	for (int i : idx) {
		for (int j = 0; j < 12; ++j) {
			A2(i, j) = gradFEy_ref2(i, j);
		}
	}
	return A2;
}

double gradFEx_ref2(const int& i, const int& j) {
	double x = 0.;
	double y = 0.;
	double grx = 0.;
	if (i == 0 || i == 4) { x = 1. / 6.;   y = 1. / 6.; }
	if (i == 1 || i == 5) { x = 2. / 3.;   y = 1. / 6.; }
	if (i == 2 || i == 6) { x = 1. / 6.;   y = 2. / 3.; }
	if (i == 3 || i == 7) { x = 1. / 3.;   y = 1. / 3.; }

	if (j == 0 && i < 4 || j == 6 && i > 3) { grx = 4 * x + 4 * y - 3; }
	if (j == 1 && i < 4 || j == 7 && i > 3) { grx = 4 * x - 1; }
	if (j == 2 && i < 4 || j == 8 && i > 3) { grx = 0; }
	if (j == 3 && i < 4 || j == 9 && i > 3) { grx = 4 - 8 * x - 4 * y; }
	if (j == 4 && i < 4 || j == 10 && i > 3) { grx = 4 * y; }
	if (j == 5 && i < 4 || j == 11 && i > 3) { grx = -4 * y; }
	return grx;
}

double gradFEy_ref2(const int& i, const int& j) {
	double x = 0.;
	double y = 0.;
	double gry = 0.;
	if (i == 0 || i == 4) { x = 1. / 6.;   y = 1. / 6.; }
	if (i == 1 || i == 5) { x = 2. / 3.;   y = 1. / 6.; }
	if (i == 2 || i == 6) { x = 1. / 6.;   y = 2. / 3.; }
	if (i == 3 || i == 7) { x = 1. / 3.;   y = 1. / 3.; }

	if (j == 0 && i < 4 || j == 6 && i > 3) { gry = 4 * x + 4 * y - 3; }
	if (j == 1 && i < 4 || j == 7 && i > 3) { gry = 0; }
	if (j == 2 && i < 4 || j == 8 && i > 3) { gry = 4 * y - 1; }
	if (j == 3 && i < 4 || j == 9 && i > 3) { gry = -4 * x; }
	if (j == 4 && i < 4 || j == 10 && i > 3) { gry = 4 * x; }
	if (j == 5 && i < 4 || j == 11 && i > 3) { gry = 4 - 4 * x - 8 * y; }
	return gry;
}

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
//double basis_2(const int& i, const int& j) {
//	double x = 0;
//	double y = 0;
//	double bas = 0;
//	if (i == 0 || i == 3) { x = 0.5; y = 0.; }
//	if (i == 1 || i == 4) { x = 0.5;   y = 0.5; }
//	if (i == 2 || i == 5) { x = 0.;   y = 0.5; }
//
//	if (j == 0 && i < 3 || j == 6 && i > 2) { bas = (1 - x - y) * (1 - 2 * x - 2 * y); }
//	if (j == 1 && i < 3 || j == 7 && i > 2) { bas = x * (2 * x - 1); }
//	if (j == 2 && i < 3 || j == 8 && i > 2) { bas = y * (2 * y - 1); }
//	if (j == 3 && i < 3 || j == 9 && i > 2) { bas = 4 * x * (1 - x - y); }
//	if (j == 4 && i < 3 || j == 10 && i > 2) { bas = 4 * x * y; }
//	if (j == 5 && i < 3 || j == 11 && i > 2) { bas = 4 * y * (1 - x - y); }
//	return bas;
//}
//
//mat basis_Assem(vector<int>& idx) {
//	mat A0(6, 12);
//
//	for (int i : idx) {
//		for (int j = 0; j < 12; ++j) {
//			A0(i, j) = basis_2(i, j);
//		}
//	}
//	return A0;
//}
//
//mat gradFEx_ref_Assem(vector<int>& idx) {
//	mat A1(6, 12);
//
//	for (int i : idx) {
//		for (int j = 0; j < 12; ++j) {
//			A1(i, j) = gradFEx_ref2(i, j);
//		}
//	}
//	return A1;
//}
//mat gradFEy_ref_Assem(vector<int>& idx) {
//	mat A2(6, 12);
//
//	for (int i : idx) {
//		for (int j = 0; j < 12; ++j) {
//			A2(i, j) = gradFEy_ref2(i, j);
//		}
//	}
//	return A2;
//}
//
//double gradFEx_ref2(const int& i, const int& j) {
//	double x = 0;
//	double y = 0;
//	double grx = 0;
//	if (i == 0 || i == 3) { x = 0.5; y = 0.; }
//	if (i == 1 || i == 4) { x = 0.5;   y = 0.5; }
//	if (i == 2 || i == 5) { x = 0.;   y = 0.5; }
//
//	if (j == 0 && i < 3 || j == 6 && i > 2) { grx = 4 * x + 4 * y - 3; }
//	if (j == 1 && i < 3 || j == 7 && i > 2) { grx = 4 * x - 1; }
//	if (j == 2 && i < 3 || j == 8 && i > 2) { grx = 0; }
//	if (j == 3 && i < 3 || j == 9 && i > 2) { grx = 4 - 8 * x - 4 * y; }
//	if (j == 4 && i < 3 || j == 10 && i > 2) { grx = 4 * y; }
//	if (j == 5 && i < 3 || j == 11 && i > 2) { grx = -4 * y; }
//	return grx;
//}
//
//double gradFEy_ref2(const int& i, const int& j) {
//	double x = 0;
//	double y = 0;
//	double gry = 0;
//	if (i == 0 || i == 3) { x = 0.5; y = 0.; }
//	if (i == 1 || i == 4) { x = 0.5;   y = 0.5; }
//	if (i == 2 || i == 5) { x = 0.;   y = 0.5; }
//
//	if (j == 0 && i < 3 || j == 6 && i > 2) { gry = 4 * x + 4 * y - 3; }
//	if (j == 1 && i < 3 || j == 7 && i > 2) { gry = 0; }
//	if (j == 2 && i < 3 || j == 8 && i > 2) { gry = 4 * y - 1; }
//	if (j == 3 && i < 3 || j == 9 && i > 2) { gry = -4 * x; }
//	if (j == 4 && i < 3 || j == 10 && i > 2) { gry = 4 * x; }
//	if (j == 5 && i < 3 || j == 11 && i > 2) { gry = 4 - 4 * x - 8 * y; }
//	return gry;
//}