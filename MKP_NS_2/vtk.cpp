#include <iostream>
#include <fstream>
#include "vtk.h"
#include <iomanip> // pro std::setprecision

void VTKfile(const std::string& filename, const string& meshFile, MatrixVectorResult& meshResult, Eigen::VectorXd& solution_x) {
    MatrixVectorResult result = meshResult; //naètení sítì
    vector<vector<double>> points = result.matrix3;
    vector<vector<int>> triangles = result.matrix1;
    int nTri = result.integerResult;
    Eigen::VectorXd x = solution_x;

    // Otevøení souboru pro zápis
    ofstream vtkFile(filename.c_str());

    // Zápis hlavièky VTK souboru
    vtkFile << "# vtk DataFile Version 3.0" << endl;
    vtkFile << "    " << endl;
    vtkFile << "ASCII" << endl;
    vtkFile << "DATASET UNSTRUCTURED_GRID" << endl;


    int numPoints = result.matrix3.size();    // Poèet bodù

    vtkFile << "POINTS " << numPoints << " float" << endl;  // Zápis informací o bodech
    //vtkFile << 0 << " " << 0 << " " << 0 << endl;
    for (int i = 0; i < (numPoints); i++) {
        //cout << points[i][0] << " " << points[i][1] << " " << points[i][2] << endl;
        vtkFile << points[i][0] << " " << points[i][1] << " " << points[i][2] << endl;
    }

    // Zápis informací o buòkách 
    vtkFile << "CELLS " << nTri << " " << nTri * 7 << endl;         // ZMENA /////////////
    //vtkFile << numPoints;
    for (int i = 0; i < nTri; ++i) {
        vtkFile << 6;
        for (int j = 0; j < 6; ++j) {
            vtkFile << " " << triangles[i][j] - 1;
        }
        vtkFile << endl;
    }


    // Zápis typu buòky
    vtkFile << "CELL_TYPES " << nTri << endl;
    for (int i = 0; i < nTri; ++i) {
        vtkFile << 22 << endl;
    }
    // Zápis hodnot øešení 
    vtkFile << "POINT_DATA " << numPoints << endl;

    //vtkFile << "SCALARS scalars float 1" << endl;
    //vtkFile << "LOOKUP_TABLE default" << endl;
    //for (int i = 0; i < numPoints; ++i) {
    //    vtkFile << x[i] << endl;
    //}

    // Nastavit počet desetinných míst
    vtkFile << std::fixed << std::setprecision(15);


    vtkFile << "VECTORS velocity double" << endl;
    for (int i = 0; i < numPoints; ++i) {
        vtkFile << x[i] << " " << x[i + numPoints] << " " << 0.0 << endl;
    }


    // Závìr souboru
    vtkFile.close();

    cout << "Soubor " << filename << " byl vytvořen." << endl;
}