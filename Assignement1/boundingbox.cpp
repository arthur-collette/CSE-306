#pragma once

#include <vector>
#include "geometry.cpp"

using namespace std;

void find_interesction(Vector N, Ray ray, vector<double> &min_interval, vector<double> &max_interval, Vector min, Vector max){
    
    double t1 = dot(min-ray.O, N) / dot(ray.u, N);
    double t2 = dot(max-ray.O, N) / dot(ray.u, N);
    if (t1 < t2){
        min_interval.push_back(t1);
        max_interval.push_back(t2);
    }
    else{
        min_interval.push_back(t2);
        max_interval.push_back(t1);
    }
}

class BoundingBox{
    public:
        Vector max;
        Vector min;
        
        BoundingBox(){}

        bool intersect(Ray &ray, double& distance){
            vector<double> min_interval;
            vector<double> max_interval;

            find_interesction(Vector(1,0,0), ray, min_interval, max_interval, min, max);
            find_interesction(Vector(0,1,0), ray, min_interval, max_interval, min, max);
            find_interesction(Vector(0,0,1), ray, min_interval, max_interval, min, max);

            if (min_interval.size() == 0){return false;}
            
            auto biggest = min_element(begin(max_interval), end(max_interval));
            auto smallest = max_element(begin(min_interval), end(min_interval));

            if((*biggest) > (*smallest) && (*smallest) > 0){
                distance = (*smallest);
                return true;
                }
            return false;
        }
    private:
};


