#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "PolygonalMesh.hpp"       
#include "Utils.hpp"               
#include "UCDUtilities.hpp"        
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
using namespace TriangularLibrary; 

// Funzione che calcola la distanza euclidea tra due punti 2D.
double computeDistance2D(const Vector2d& p1, const Vector2d& p2) {
    return (p2 - p1).norm();
}

// Funzione che calcola l'area di un poligono (in 2D) utilizzando la formula di Gauss.
double computePolygonArea(const vector<Vector2d>& vertices) {
    double area = 0.0;
    int n = vertices.size();
    for (int i = 0; i < n; i++) {
        const Vector2d& p1 = vertices[i];
        const Vector2d& p2 = vertices[(i + 1) % n];
        area += p1.x() * p2.y() - p2.x() * p1.y();
    }
    return fabs(area) / 2.0;
}

// Funzione che converte la matrice di punti in un vettore di punti 2D
vector<Vector2d> get2DPoints(const MatrixXd& points) {
    vector<Vector2d> pts;
    int numPoints = points.cols();
    pts.reserve(numPoints);
    for (int i = 0; i < numPoints; ++i) {
        double x = points(0, i);
        double y = points(1, i);
        pts.push_back(Vector2d(x, y));
    }
    return pts;
}

int main() {
    PolygonalMesh mesh;
    if (!ImportMesh(mesh)) {
        cerr << "File not found or error importing mesh" << endl;
        return 1;
    }

    const double epsilon = 1e-9;

    for (const auto& entry : mesh.MarkerCell0Ds) {
        int marker = entry.first;
        // Il marker deve essere non negativo.
        assert(marker >= 0);
    }

    for (const auto& entry : mesh.MarkerCell1Ds) {
        int marker = entry.first;
        assert(marker >= 0);
    }

    MatrixXd points = mesh.Cell0DsCoordinates;  
    vector<Vector2d> points2D = get2DPoints(points);  // Prendo solo le coordinate X e Y

    int numEdges = mesh.Cell1DsExtrema.cols();
    for (int e = 0; e < numEdges; ++e) {
        int idxOrigin = mesh.Cell1DsExtrema(0, e);
        int idxEnd    = mesh.Cell1DsExtrema(1, e);
        // Recupero dei punti 2D corrispondenti agli indici
        Vector2d pOrigin = points2D[idxOrigin];
        Vector2d pEnd    = points2D[idxEnd];
        double length = computeDistance2D(pOrigin, pEnd);
        assert(length > epsilon);
    }
    for (const auto& polygon : mesh.Cell2Ds) {
        vector<Vector2d> polyVertices;
        for (unsigned int idx : polygon.Vertices) {
            polyVertices.push_back(points2D[idx]);
        }
        double area = computePolygonArea(polyVertices);
        assert(area > epsilon);
    }

    cout << "Tutti i test sul mesh sono stati superati correttamente." << endl;

    Gedim::UCDUtilities utilities;

    vector<Gedim::UCDProperty<double>> pointProperties(1);
    pointProperties[0].Label = "Marker";
    pointProperties[0].UnitLabel = "-";
    pointProperties[0].NumComponents = 1;

    vector<double> pointMarkers(mesh.NumCell0Ds, 0.0);
    for (const auto& entry : mesh.MarkerCell0Ds) {
        int marker = entry.first;
        for (unsigned int idx : entry.second) {
            if (idx < pointMarkers.size())
                pointMarkers[idx] = marker;
        }
    }
    pointProperties[0].Data = pointMarkers.data();

    utilities.ExportPoints("./Cell0Ds.inp", mesh.Cell0DsCoordinates, pointProperties);

    vector<Gedim::UCDProperty<double>> segProperties(1);
    segProperties[0].Label = "Marker";
    segProperties[0].UnitLabel = "-";
    segProperties[0].NumComponents = 1;

    vector<double> segMarkers(mesh.NumCell1Ds, 0.0);
    for (const auto& entry : mesh.MarkerCell1Ds) {
        int marker = entry.first;
        for (unsigned int idx : entry.second) {
            if (idx < segMarkers.size())
                segMarkers[idx] = marker;
        }
    }
    segProperties[0].Data = segMarkers.data();

    utilities.ExportSegments("./Cell1Ds.inp", mesh.Cell0DsCoordinates, mesh.Cell1DsExtrema, {}, segProperties);

    cout << "Esportazione completata." << endl;

    return 0;
}

