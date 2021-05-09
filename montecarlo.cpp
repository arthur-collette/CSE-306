#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>

static std::default_random_engine engine(10) ; // random seed = 10 
static std::uniform_real_distribution<double> uniform(0, 1);

void boxMuller(double stdev , double &x, double &y) { 
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    x = sqrt(-2. * log(r1))*cos(2.*M_PI*r2)*stdev; 
    y = sqrt(-2. * log(r1))*sin(2.*M_PI*r2)*stdev;
}

double gauss(double& x, double& stdev){
    return (1/(stdev * sqrt(2*M_PI))) * exp(-pow(x,2.) / (2*pow(stdev,2.)));
}

int main(){
    int N = 10000;
    double stdev = 1.;
    double result = 0.;
    for (int i=0; i<N; i++){
        double w,x,y,z;
        boxMuller(stdev,w,x);
        boxMuller(stdev,y,z);
        if (w < -M_PI/2 | w > M_PI/2 | x < -M_PI/2 | x > M_PI/2 || y < -M_PI/2 || y > M_PI/2){ continue; } 
        result += cos(x*y*z) / (gauss(w,stdev)*gauss(x,stdev)*gauss(y,stdev));               
    }
    result /= N;
    std::cout<<result<<std::endl;
    return 0;
}