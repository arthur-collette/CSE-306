#pragma once
#include <iostream>
#include <vector>

#include "polygon.cpp"
#include "svg.cpp"

using namespace std;

Vector intersect(Vector &A, Vector &B, Edge &edge) {
    Vector u = edge.vertices[0];
    Vector v = edge.vertices[1];
    Vector N = normalize(Vector(v[1] - u[1], u[0] - v[0]));
    double t = dot(u - A, N) / dot(B - A, N);
    return A + (B - A)*t;
}

bool inside(Vector &P, Edge &edge) {
    Vector u = edge.vertices[0];
    Vector v = edge.vertices[1];
    Vector N = Vector(v[1] - u[1], u[0] - v[0]);
    if (dot(P - u, N) <= 0) return true;
    return false;
}

// Sutherland-Hodgman Polygon Clipping
Polygon clip(Polygon subjectPolygon, Polygon clipPolygon) {
    Polygon outPolygon;
    for (int i = 0; i < clipPolygon.vertices.size(); i++) {
        Vector u = clipPolygon.vertices[i];
        Vector v = clipPolygon.vertices[(i>0)?(i-1):(clipPolygon.vertices.size()-1)];
        Edge edge(u, v);
        outPolygon = Polygon();
        for (int j = 0; j < subjectPolygon.vertices.size(); j++) {
            Vector curVertex = subjectPolygon.vertices[j];
            Vector prevVertex = subjectPolygon.vertices[(j>0)?(j-1):(subjectPolygon.vertices.size() - 1)];
            Vector intersection = intersect(prevVertex, curVertex, edge);
            if (inside(curVertex, edge)) {
                if (!inside(prevVertex, edge)) {
                    outPolygon.vertices.push_back(intersection);
                }
                outPolygon.vertices.push_back(curVertex);
            } else if (inside(prevVertex, edge)) {
                outPolygon.vertices.push_back(intersection);
            }
        }
        subjectPolygon = outPolygon;
    }
    return outPolygon;
}