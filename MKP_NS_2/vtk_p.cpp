#include <iostream>
#include <fstream>
#include "vtk_p.h"
#include <iomanip> // pro std::setprecision

void VTKfile_p(const std::string& filename, const string& meshFile, MatrixVectorResult& meshResult, Eigen::VectorXd& solution_x) {
    MatrixVectorResult result = meshResult; //naètení sítì
    vector<vector<double>> points = result.matrix3;
    vector<vector<int>> triangles = result.matrix1;
    int nTri = result.integerResult;
    Eigen::VectorXd x = solution_x;

    vector<int> M(points.size(), 0);
    vector<int> MM(points.size(), 0);
    int kk = 0;
    for (int k = 0; k < nTri; k++) {
        int A = result.matrix1[k][0] - 1;					// Oznaèení vrcholù trojúhelníku k
        int B = result.matrix1[k][1] - 1;
        int C = result.matrix1[k][2] - 1;
        vector<int> verTri = { A, B, C };
        for (int i = 0; i < 3; i++) {
            M[verTri[i]] = 1;
        }

    }
    for (int k = 0; k < points.size(); k++) {
        if (M[k] == 1) {
            MM[k] = kk;
            kk++;
        }
    }
    //cout << "M" << endl;
    //for (int i = 0; i < points.size(); ++i) {
    //    cout << M[i] << endl;
    //}

    //cout << "MM" << endl;
    //for (int i = 0; i < points.size(); ++i) {
    //    cout << MM[i] << endl;
    //}
    vector<vector<double>> newPoints;

    for (int i = 0; i < points.size(); ++i) {
        if (M[i] == 1) {
            newPoints.push_back(points[i]);
        }
    }

    //// Výpis hodnot
    //std::cout << "newPoints obsahuje:\n";
    //for (const auto& row : newPoints) {
    //    for (const auto& val : row) {
    //        std::cout << val << " ";
    //    }
    //    std::cout << std::endl;
    //}

    vector<vector<int>> newTriangles(nTri, std::vector<int>(3, 0));
    for (int k = 0; k < nTri; k++) {
        for (int i = 0; i < 3; i++) {
            newTriangles[k][i] = MM[triangles[k][i] - 1];
            /*cout << "k " << k << " i  " << i << "  triangles[k][i] " << triangles[k][i] << "  MM[triangles[k][i]]  " << MM[triangles[k][i]] << endl;*/
        }
    }
    //// Výpis hodnot
    //std::cout << "newtriangles obsahuje:\n";
    //for (const auto& row : newTriangles) {
    //    for (const auto& val : row) {
    //        std::cout << val << " ";
    //    }
    //    std::cout << std::endl;
    //}

    // Otevøení souboru pro zápis
    ofstream vtkFile(filename.c_str());

    // Zápis hlavièky VTK souboru
    vtkFile << "# vtk DataFile Version 3.0" << endl;
    vtkFile << "    " << endl;
    vtkFile << "ASCII" << endl;
    vtkFile << "DATASET UNSTRUCTURED_GRID" << endl;


    int numPoints = newPoints.size();    // Poèet bodù

    vtkFile << "POINTS " << numPoints << " float" << endl;  // Zápis informací o bodech
    //vtkFile << 0 << " " << 0 << " " << 0 << endl;
    for (int i = 0; i < (numPoints); i++) {
        //cout << points[i][0] << " " << points[i][1] << " " << points[i][2] << endl;
        vtkFile << newPoints[i][0] << " " << newPoints[i][1] << " " << newPoints[i][2] << endl;
    }

    // Zápis informací o buòkách 
    vtkFile << "CELLS " << nTri << " " << nTri * 4 << endl;         // ZMENA /////////////
    //vtkFile << numPoints;
    for (int i = 0; i < nTri; ++i) {
        vtkFile << 3;
        for (int j = 0; j < 3; ++j) {
            vtkFile << " " << newTriangles[i][j];
        }
        vtkFile << endl;
    }


    // Zápis typu buòky
    vtkFile << "CELL_TYPES " << nTri << endl;
    for (int i = 0; i < nTri; ++i) {
        vtkFile << 5 << endl;
    }
    // Zápis hodnot øešení 
    vtkFile << "POINT_DATA " << x.size() << endl;

    // Nastavit počet desetinných míst
    vtkFile << std::fixed << std::setprecision(15);

    vtkFile << "SCALARS scalars float 1" << endl;
    vtkFile << "LOOKUP_TABLE default" << endl;
    for (int i = 0; i < x.size(); ++i) {
        vtkFile << x[i] << endl;
    }



    //vtkFile << "VECTORS displacement double" << endl;
    //for (int i = 0; i < numPoints; ++i) {
    //    vtkFile << x[i] << " " << x[i + numPoints] << " " << 0.0 << endl;
    //}


    // Závìr souboru
    vtkFile.close();

    std::cout << "Soubor " << filename << " byl vytvořen." << endl;
}