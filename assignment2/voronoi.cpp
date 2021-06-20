#pragma once
#include <iostream>
#include <vector>

#include "polygon.cpp"
#include "svg.cpp"

Vector intersect_V(Vector A, Vector B, Vector Pi, Vector Pj){
    Vector N = (Pi + Pj)/2;
    double t = dot(N - A, Pi - Pj)/dot(B - A, Pi - Pj);
    return A + (B - A) * t;
}

bool inside_V(Vector P, Vector Pi, Vector Pj){
    Vector N = (Pi + Pj)/2;
    if (dot(P - N, Pj - Pi) < 0) return true;
    return false;
}

vector<Polygon> voronoi_diagram(Polygon points, Polygon box){
    std::vector<Polygon> diagram;
    Polygon cell;
    for(int i = 0; i < points.vertices.size(); i++) {
        Vector Pi = points.vertices[i];
        Polygon cell = box;
        for(int j = 0; j < points.vertices.size(); j++) {
            if(i!=j){
                Polygon outPolygon = Polygon();
                Vector Pj = points.vertices[j];
                for(int k = 0; k < cell.vertices.size(); k++) {
                    Vector curVertex = cell.vertices[k];
                    Vector prevVertex = cell.vertices[(k>0)?(k-1):(cell.vertices.size()-1)];
                    Vector intersection = intersect_V(prevVertex, curVertex, Pi, Pj);
                    if (inside_V(curVertex, Pi, Pj)) {
                        if (!inside_V(prevVertex, Pi, Pj)) {
                            outPolygon.vertices.push_back(intersection);
                        }
                        outPolygon.vertices.push_back(curVertex);
                    } else if (inside_V(prevVertex, Pi, Pj)) {
                        outPolygon.vertices.push_back(intersection);
                    }
                }
                cell = outPolygon;
            } 
        }
        diagram.push_back(cell);
    }
    return diagram;
}


Vector intersect_PV(Vector A, Vector B, Vector Pi, Vector Pj, double weight_i, double weight_j){
    Vector N = (Pi + Pj)/2 + (Pj - Pi) * (weight_i - weight_j)/(2*pow(norm(Pi - Pj), 2));
    double t = dot(N - A, Pi - Pj)/dot(B - A, Pi - Pj);
    return A + (B - A) * t;
}

bool inside_PV(Vector P, Vector Pi, Vector Pj, double weight_i, double weight_j){
    Vector N = (Pi + Pj)/2 + (Pj - Pi) * (weight_i - weight_j)/(2*pow(norm(Pi - Pj), 2));
    if (dot(P - N, Pj - Pi) < 0) return true;
    return false;
}

vector<Polygon> voronoi_diagram(Polygon points, Polygon box, vector<double> weights){
    std::vector<Polygon> diagram;
    Polygon cell;
    for(int i = 0; i < points.vertices.size(); i++) {
        Vector Pi = points.vertices[i];
        double weight_i = weights[i];
        Polygon cell = box;
        for(int j = 0; j < points.vertices.size(); j++) {
            if(i!=j){
                Polygon outPolygon = Polygon();
                Vector Pj = points.vertices[j];
                double weight_j = weights[j];
                for(int k = 0; k < cell.vertices.size(); k++) {
                    Vector curVertex = cell.vertices[k];
                    Vector prevVertex = cell.vertices[(k>0)?(k-1):(cell.vertices.size()-1)];
                    Vector intersection = intersect_PV(prevVertex, curVertex, Pi, Pj, weight_i, weight_j);
                    if (inside_PV(curVertex, Pi, Pj, weight_i, weight_j)) {
                        if (!inside_PV(prevVertex, Pi, Pj, weight_i, weight_j)) {
                            outPolygon.vertices.push_back(intersection);
                        }
                        outPolygon.vertices.push_back(curVertex);
                    } else if (inside_PV(prevVertex, Pi, Pj, weight_i, weight_j)) {
                        outPolygon.vertices.push_back(intersection);
                    }
                }
                cell = outPolygon;
            } 
        }
        diagram.push_back(cell);
    }
    return diagram;
}