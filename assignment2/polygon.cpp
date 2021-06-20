#pragma once

#include <vector>

#include "vector.cpp"

using namespace std;

struct Edge{
    vector<Vector> vertices;

    Edge(){}
    Edge(Vector v1, Vector v2){
        vertices.push_back(v1);
        vertices.push_back(v2);
    }
};

class Polygon{
    public: 
        vector<Vector> vertices;

        Polygon(){};
        Polygon(vector<Vector> vertices){
            this->vertices = vertices;
        }
    
    private:
};