
#include <iostream>
#include "PolygonalMesh.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace PolygonalLibrary;


bool TestEdges(const PolygonalMesh &mesh) {
   
	const double epsilon = 1e-10;

    for (unsigned int cell_id = 0; cell_id < mesh.NumCell1Ds; ++cell_id) {
        unsigned int originId = mesh.Cell1DsExtrema(0, cell_id);
        unsigned int endId = mesh.Cell1DsExtrema(1, cell_id);
		
		Vector3d origin = mesh.Cell0DsCoordinates.col(originId);
        Vector3d end = mesh.Cell0DsCoordinates.col(endId);

        double length = (end - origin).norm();
        if (length < epsilon) {
            cerr << "Errore: il segmento " << cell_id << " ha lunghezza nulla" << endl;
            return false;
        }
    }
    return true; 
}



bool TestArea(const PolygonalMesh &mesh) {
    
	const double epsilon = 1e-10; 

    for (unsigned int cell_id = 0; cell_id < mesh.NumCell2Ds; ++cell_id) {
        const auto& vertices = mesh.Cell2DsVertices[cell_id]; // Poligono corrente
		if (vertices.size() < 3) continue;

        // Calcolo l'area del poligono
        double area = 0.0;
		int n = vertices.size();
		for (int i = 0; i < n; ++i) {
			unsigned int k = (i + 1) % n;
			const auto& v1 = mesh.Cell0DsCoordinates.col(vertices[i]);
            const auto& v2 = mesh.Cell0DsCoordinates.col(vertices[k]);
			area += v1.x() * v2.y() - v1.y() * v2.x();
    }
    return std::abs(area) / 2.0;
        
    if (area < epsilon) {
            cerr << "Errore: poligono con indice " << cell_id << " ha area nulla" << endl;
            return false;
        }
    }
    return true;
}


int main()
{
    PolygonalMesh mesh;

    if(!ImportMesh(mesh))
    {
        cerr << "file not found" << endl;
        return 1;
    }

    Gedim::UCDUtilities utilities;
    {
        vector<Gedim::UCDProperty<double>> cell0Ds_properties(1);

        cell0Ds_properties[0].Label = "Marker";
        cell0Ds_properties[0].UnitLabel = "-";
        cell0Ds_properties[0].NumComponents = 1;

        vector<double> cell0Ds_marker(mesh.NumCell0Ds, 0.0);
        for(const auto &m : mesh.MarkerCell0Ds)
            for(const unsigned int id: m.second)
                cell0Ds_marker.at(id) = m.first;

        cell0Ds_properties[0].Data = cell0Ds_marker.data();

        utilities.ExportPoints("./Cell0Ds.inp",
                               mesh.Cell0DsCoordinates,
                               cell0Ds_properties);
    }

    {

        vector<Gedim::UCDProperty<double>> cell1Ds_properties(1);

        cell1Ds_properties[0].Label = "Marker";
        cell1Ds_properties[0].UnitLabel = "-";
        cell1Ds_properties[0].NumComponents = 1;

        vector<double> cell1Ds_marker(mesh.NumCell1Ds, 0.0);
        for(const auto &m : mesh.MarkerCell1Ds)
            for(const unsigned int id: m.second)
                cell1Ds_marker.at(id) = m.first;

        cell1Ds_properties[0].Data = cell1Ds_marker.data();

        utilities.ExportSegments("./Cell1Ds.inp",
                                 mesh.Cell0DsCoordinates,
                                 mesh.Cell1DsExtrema,
                                 {},
                                 cell1Ds_properties);
    }
	
	
	/*
	{
		vector<Gedim::UCDProperty<double>> cell0Ds_properties(1);

        cell0Ds_properties[0].Label = "Marker";
        cell0Ds_properties[0].UnitLabel = "-";
        cell0Ds_properties[0].NumComponents = 1;

        vector<double> cell0Ds_marker(mesh.NumCell0Ds, 0.0);
        for(const auto &m : mesh.MarkerCell0Ds)
            for(const unsigned int id: m.second)
                cell0Ds_marker.at(id) = m.first;

        cell0Ds_properties[0].Data = cell0Ds_marker.data();
		
		vector<Gedim::UCDProperty<double>> cell2Ds_properties(1);

		cell2Ds_properties[0].Label = "Marker";
		cell2Ds_properties[0].UnitLabel = "-";
		cell2Ds_properties[0].NumComponents = 1;

		vector<double> cell2Ds_marker(mesh.NumCell2Ds, 0.0);
		for(const auto &m : mesh.MarkerCell2Ds)
			for(const unsigned int id: m.second)
				cell2Ds_marker.at(id) = m.first;

		cell2Ds_properties[0].Data = cell2Ds_marker.data();

		utilities.ExportPolygons("./Cell2Ds.inp",
								 mesh.Cell0DsCoordinates,
								 mesh.Cell2DsVertices,
								 cell0Ds_properties,
								 cell2Ds_properties);
 
    }*/
	
	// marker control
    map<unsigned int, list<unsigned int>> MarkerCell0Ds_true = {
    {1, {0}},
    {2, {1}},
    {3, {2}},
    {4, {3}},
    {5, {6, 16, 24}},
    {6, {7, 17, 22, 78}},
    {7, {8, 20, 23, 52, 59}},
    {8, {5, 15, 21, 26, 92}}
    };
    
    if ( mesh.MarkerCell0Ds == MarkerCell0Ds_true) {
        cout << "Marker celle 0D corretti" << endl;
    } else {
        cout << "Marker celle 0D errati" << endl;
    }
    
    map<unsigned int, list<unsigned int>> MarkerCell1Ds_true = {
	{5, {8,19,22,28}},
    {6, {6, 23, 26,126,127}},
    {7, {14,17,24,79,92, 93}},
    {8, {11,25,29,30,159,160}},
    };
    
    if ( mesh.MarkerCell1Ds == MarkerCell1Ds_true) {
        cout << "Marker celle 1D corretti" << endl;
    } else {
        cout << "Marker celle 1D errati" << endl;
    }
    
    map<unsigned int, list<unsigned int>> MarkerCell2Ds_true = {};
    
    if ( mesh.MarkerCell2Ds == MarkerCell2Ds_true) {
        cout << "Marker celle 2D corretti" << endl;
    } else {
        cout << "Marker celle 2D errati" << endl;
    }
	
	
	// Test sui lati del poligoni
    bool result_edges = TestEdges(mesh);
    if (result_edges) {
        cout << "Tutti i segmenti hanno lunghezza non nulla" << endl;
    }
    
    // controllo sulle aree
    bool result_area = TestArea(mesh);
    if (result_area) {
        cout << "Tutti i poligoni hanno area non nulla" << endl;
    }

	
    return 0;
}