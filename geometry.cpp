#pragma once
#include "ray.cpp"
#include <iostream>  

using namespace std;

struct Intersection
{
    bool result;
    Vector P;
    Vector N;
    double t;
    Vector albedo;

    Intersection(){
        result = false;
    }
    Intersection(bool b){
        result = b;
    }
    Intersection(bool result, Vector P, Vector N, double t, Vector albedo){
        this->result = result;
        this->P = P;
        this->N = N;
        this->t = t;
        this->albedo = albedo;
    }
};

class Geometry{
    public: 

        virtual Intersection intersect(Ray ray) = 0; 

        Vector get_albedo(){return albedo;}
        bool is_mirror(){return mirror;}

        Vector albedo;
        bool mirror;
        bool refract;
        double n;

    private:
        
};