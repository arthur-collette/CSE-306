#pragma once
#include "vector.cpp"
#include <iostream>

class Ray{
    public:
        Vector O;
        Vector u;

        Ray(Vector O, Vector u){
            this->O = O;
            this->u = u;
        }

    private:

};