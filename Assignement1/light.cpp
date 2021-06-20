#pragma once
#include <iostream>
#include "vector.cpp"

class Light{

    public:
        Vector S;
        double I;

        Light(){
            this->S = Vector(0,0,0);
            this->I = 0;
        }

        Light(Vector S, double I){
            this->S = S;
            this->I = I;
        }
};
