#pragma once
#include "ray.cpp"

class Camera{
    public:
        Vector Q;
        int H;
        int W;
        double alpha;

        Camera(Vector Q, int H, int W, double alpha){
            this->Q = Q;
            this->H = H;
            this->W = W;
            this->alpha= alpha;
        }
};

