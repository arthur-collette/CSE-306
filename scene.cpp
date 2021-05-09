#pragma once
#include "sphere.cpp"
#include "light.cpp"

#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>
#include <math.h>
#include <vector>

static default_random_engine engine(10);
static uniform_real_distribution<double> uniform(0, 1);

using namespace std;

int V(Ray ray, Vector w, double d, vector<Geometry *> arrayG) {
    double epsilon = 0.01;
    for(vector<Geometry*>::iterator it = arrayG.begin() ; it != arrayG.end(); ++it){
        Intersection in = (*it)->intersect(ray);
        if (in.result){
            if (in.t < d) return 0.;
        }
    }
    return 1;
}

Vector random_cos(const Vector &N) {

    double r1 = uniform(engine);
    double r2 = uniform(engine);

    double x = cos(2 * M_PI * r1) * sqrt(1 - r2);
    double y = sin(2 * M_PI * r1) * sqrt(1 - r2);
    double z = sqrt(r2);

    Vector T1;
    double min_component = N[0];
    double min_index = 0;
    for (int i=1;i<3;i++){
        if (N[i] < min_component){ 
            min_component = N[i];
            min_index = i;
        }
    }
    if (min_index == 0){ T1 = Vector(0.,N[2],-N[1]); }
    else if (min_index == 1){ T1 = Vector(N[2],0.,-N[0]); }
    else { T1 = Vector(N[1],-N[0],0.); }

    T1 = divide(T1,norm(T1));

    Vector T2 = cross(N, T1);

    return divide(T1 * x + T2 * y + N * z, norm(T1 * x + T2 * y + N * z));
}


class Scene{
    public:
        Scene(vector<Geometry *> arrayG, Light light){
            this->arrayG = arrayG;
            this->light = light;
        }       

        Vector getColor(const Ray& ray, int ray_depth){
            double dist = 0;
            Intersection best;
            double geo_index;
            double epsilon = 0.01;

            best = findclostest(ray, geo_index);
            

            if (ray_depth < 0){
                return Vector(0., 0., 0.);
            }

            if (arrayG[geo_index]->mirror){
                Ray relfected = Ray(best.P + best.N*epsilon, ray.u - best.N * 2*dot(divide(ray.u,norm(ray.u)),best.N)); 
                return getColor(relfected , ray_depth - 1);
            }
            if (arrayG[geo_index]->refract){
                
                double n1;
                double n2;
				Vector N_refra;

                if (dot(best.N,ray.u) > 0.) {//exiting
					N_refra = best.N * (-1.);
                    n1 = arrayG[geo_index]->n;
                    n2 = 1.;
				}
                else {//entering
				    N_refra = best.N;
                    n1 = 1.;
                    n2 = arrayG[geo_index]->n;
				}
                
                double ratio = n1/n2;
                double uN = dot(ray.u, N_refra);
                double x = 1 - ratio*ratio * (1 - pow(uN,2));   

                Vector omega_t = (ray.u - (N_refra * uN)) * ratio;
                Vector omega_n = N_refra * (-1) * sqrt(x);
                Vector omega = omega_t + omega_n;

				Ray refracted_ray = Ray(best.P - N_refra * epsilon, omega); 

				return getColor(refracted_ray, ray_depth - 1);

            }
            else{

                //direct lightening
                double d = norm(light.S - best.P);
                Vector w = divide(light.S - best.P, d);
                Ray new_ray = Ray(best.P + best.N * epsilon, w);
                double v = V(new_ray, w, d, arrayG);

                Vector color = best.albedo * ((light.I * v * max(dot(best.N, w), 0.))/ (4*pow(d,2) * pow(M_PI,2)));
                
                //indirect lighting
                Ray random_ray = Ray(best.P + best.N * epsilon, random_cos(best.N));
                Vector indirect_contribution = best.albedo * getColor(random_ray, ray_depth - 1);
                color = (color + indirect_contribution) * 0.5;

                return color;
            }
        }

        Intersection findclostest(Ray ray, double &geo_index){
            Intersection best;
            double dist = 0;
            for(int i=0; i < arrayG.size(); i++){

                Intersection inter = arrayG[i]->intersect(ray);
                
                if (inter.result){
                    if (dist == 0)
                    {
                        dist = norm(inter.P - ray.O);
                        best = inter;
                        geo_index = i;
                    }
                    else {
                        if (norm(inter.P - ray.O) < dist){
                            dist = norm(inter.P - ray.O);
                            best = inter;
                            geo_index = i;
                        }
                    }
                }
            }
            
            return best;
        }

    private:
        vector<Geometry *> arrayG;
        Light light;
};

