#include "mesh_processor.h"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <tuple>


MatrixVectorResult GetTriangles(const string& meshFile) {
    MatrixVectorResult result;
    ifstream file(meshFile);

    if (!file.is_open()) {
        cerr << "Nelze otevrit soubor " << meshFile << endl;
        return {};
    }

    string line;
    bool isReadingElements = false;
    vector < vector<double>> nodes_coordinates;

    int i = 0;
    int numNode;

    while (getline(file, line)) {
        if (line.find("$Nodes") != string::npos) {
            isReadingElements = true;
            std::getline(file, line); // P�esko� prvn� ��dek s popisem
            istringstream iss2(line);
            iss2 >> numNode;
            std::cout << numNode << endl;
            nodes_coordinates.resize(numNode);

            continue;
        }

        if (line.find("$EndNodes") != string::npos) {
            break;
        }

        if (isReadingElements) {
            istringstream iss(line);
            int nodeId;
            double cooX, cooY, cooZ;
            vector<double> coordinates(3);
            iss >> nodeId >> cooX >> cooY >> cooZ;


            /* if (i > 0) {*/

            nodes_coordinates[i] = { cooX, cooY, cooZ };
            /*}*/
            i++;
        }
    }

    isReadingElements = false; // Resetujeme na��t�n� prvk�
    file.clear(); // Resetujeme stav souboru
    file.seekg(0, std::ios::beg); // Nastav�me pozici na za��tek souboru




    vector<int> element_type_line_V;      // p��prav vektor� pro �se�ky
    vector<int> element_id_line_V;
    vector<int> num_tags_line_V;
    vector<int> physical_tag_line_V;
    vector<int> elementary_tag_line_V;

    vector<int> element_type_tri_V;       // p��prav vektor� pro troj�heln�ky
    vector<int> element_id_tri_V;
    vector<int> num_tags_tri_V;
    vector<int> physical_tag_tri_V;
    vector<int> elementary_tag_tri_V;

    vector<vector<int>> triangle_vertex_M;
    vector<vector<int>> line_vertex_M;



    int skip_first = 0;

    while (getline(file, line)) {
        if (line.find("$Elements") != string::npos) {
            isReadingElements = true;
            continue;
        }

        if (line.find("$EndElements") != string::npos) {
            break;
        }

        if (isReadingElements) {
            istringstream iss(line);
            int elementType, elementId, numTags, physicalTag, elementaryTag, nodeA, nodeB, nodeC, nodeD, nodeE, nodeF;
            iss >> elementId >> elementType >> numTags >> physicalTag >> elementaryTag >> nodeA >> nodeB >> nodeC >> nodeD >> nodeE >> nodeF;

            //std::cout << elementId << " " << elementType << " " << numTags << " " << physicalTag << " " << elementaryTag << " " << nodeA << " " << nodeB << " " << nodeC << " " << std::endl;


            if (skip_first > 0) {

                if (elementType == 1) {                              // 1 je ozna�en� �se�ek
                    elementary_tag_line_V.push_back(elementaryTag);
                    element_type_line_V.push_back(elementType);
                    element_id_line_V.push_back(elementId);
                    num_tags_line_V.push_back(numTags);
                    physical_tag_line_V.push_back(physicalTag);
                    elementary_tag_line_V.push_back(elementaryTag);
                    line_vertex_M.push_back({ nodeA, nodeB, nodeC });
                }
                if (elementType == 2) {                              // 2 je ozna�en� troj�heln�k�
                    elementary_tag_tri_V.push_back(elementaryTag);
                    element_type_tri_V.push_back(elementType);
                    element_id_tri_V.push_back(elementId);
                    num_tags_tri_V.push_back(numTags);
                    physical_tag_tri_V.push_back(physicalTag);
                    elementary_tag_tri_V.push_back(elementaryTag);
                    triangle_vertex_M.push_back({ nodeA, nodeB, nodeC });
                }
                if (elementType == 8) {                              // 8 je ozna�en� kvadratick�ch �se�ek
                    elementary_tag_line_V.push_back(elementaryTag);
                    element_type_line_V.push_back(elementType);
                    element_id_line_V.push_back(elementId);
                    num_tags_line_V.push_back(numTags);
                    physical_tag_line_V.push_back(physicalTag);
                    elementary_tag_line_V.push_back(elementaryTag);
                    line_vertex_M.push_back({ nodeA, nodeB, nodeC });
                }
                if (elementType == 9) {                              // 9 je ozna�en� kvadratick�ch troj�heln�k�
                    elementary_tag_tri_V.push_back(elementaryTag);
                    element_type_tri_V.push_back(elementType);
                    element_id_tri_V.push_back(elementId);
                    num_tags_tri_V.push_back(numTags);
                    physical_tag_tri_V.push_back(physicalTag);
                    elementary_tag_tri_V.push_back(elementaryTag);
                    triangle_vertex_M.push_back({ nodeA, nodeB, nodeC, nodeD, nodeE, nodeF });
                }
            }
            skip_first++;

        }

    }


    result.matrix1 = triangle_vertex_M;          // ukl�d�n� v�stupu funkce
    result.matrix2 = line_vertex_M;
    result.matrix3 = nodes_coordinates;
    result.vector1 = physical_tag_line_V;
    result.vector2 = physical_tag_tri_V;

    int n_Tri = physical_tag_tri_V.size();      // po�et troj�heln�k�
    int n_Line = physical_tag_line_V.size();    // po�et �se�ek

    result.integerResult = n_Tri;
    result.integerResult2 = n_Line;


    file.close();
    return result;
}