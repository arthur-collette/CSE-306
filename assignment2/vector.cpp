#pragma once
#include <iostream>

using namespace std;

class Vector {
    public:
        explicit Vector(double x = 0., double y = 0., double z = 0.) {
            coords[0] = x;
            coords[1] = y;
            coords[2] = z;
        };
        Vector& operator+=(const Vector& b){
            coords[0] += b[0] ;
            coords[1] += b[1] ;
            coords[2] += b[2] ;
            return *this ;
        }
        const double& operator[] (int i) const { return coords[i]; }
        double& operator[] (int i) { return coords[i]; }

        void print(){
            cout << coords[0] << " " << coords[1] << " " << coords[2] << endl;
        }
    private:
        double coords[3];
};

Vector operator+(const Vector& a, const Vector& b){
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vector operator-(const Vector& a, const Vector& b){
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

// Vector operator*(const Vector& a, const Vector& b){
//     return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
// }

Vector operator*(const Vector& a, double b){
    return Vector(b * a[0], b *a[1], b *a[2]);
}

Vector operator/(const Vector& a, double t) {
  return Vector(a[0]/t, a[1]/t, a[2]/t);
}

double dot(const Vector& a, const Vector& b){
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double norm(const Vector& a){
    return sqrt(pow(a[0],2) + pow(a[1],2) + pow(a[2],2));
}

Vector divide(const Vector& a, double n){
    return Vector(a[0]/n, a[1]/n, a[2]/n);
}

Vector cross(const Vector& a, const Vector& b){
    return Vector(a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]);
}

Vector normalize(const Vector& a) {
  return a/sqrt(dot(a, a));
}
