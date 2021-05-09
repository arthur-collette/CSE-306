#pragma once
#include "geometry.cpp"

using namespace std;

class Sphere : public Geometry{
    private:
        Vector C;
        double R;

    public:
        Sphere(){
            this->C = Vector(0,0,0);
            this->R = 0;
            this->albedo = Vector(0,0,0);
            this->mirror = false;
            this->refract = false;
            n = 1.5;
        }

        Sphere(Vector C, double R, Vector albedo, bool mirror, bool refract){
            this->C = C;
            this->R = R;
            this->albedo = albedo;
            this->mirror = mirror;
            this->refract = refract;
            n = 1.5;
        }

        Intersection intersect(Ray ray){
            double delta = pow(dot(ray.u, (ray.O-C)),2)  - (pow(norm(ray.O-C),2) - pow(R,2));
            if (delta < 0){
                return Intersection(false);
            }

            else if (delta == 0){

                double t = dot(ray.u, (C-ray.O));

                if (t < 0){
                    return Intersection(false);
                }
                Vector P = ray.O + ray.u*t;
                return Intersection(true, P, divide(P-C, norm(P-C)), t, albedo);
            }
            else {
                double t2 = dot(ray.u, C-ray.O) + sqrt(delta);
                if (t2 < 0){
                    return Intersection(false);
                }
                double t1 = dot(ray.u, C-ray.O) - sqrt(delta);
                if (t1 >= 0){
                    Vector P = ray.O + ray.u*t1;
                    return Intersection(true, P, divide(P-C, norm(P-C)), t1, albedo);
                }
                else{
                    Vector P = ray.O + ray.u*t2;
                    return Intersection(true, P, divide(P-C, norm(P-C)), t2, albedo);
                }
            }
        }
};

