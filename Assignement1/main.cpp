#include <iostream>
#include <vector>
#include <math.h> 

#include "camera.cpp"
#include "scene.cpp"
#include "mesh.cpp"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

using namespace std;

void boxMuller(double stdev , double &x, double &y) { 
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    x = sqrt(-2. * log(r1))*cos(2.*M_PI*r2)*stdev; 
    y = sqrt(-2. * log(r1))*sin(2.*M_PI*r2)*stdev;
}


int main() {

    int W = 512;
    int H = 512;
    int scene_nbr = 0; //0 for 3 balls, 1 for cat and 2 for 1 ball
    std::vector<unsigned char> image(W*H * 3, 0);

    TriangleMesh triangleMesh;
    triangleMesh.readOBJ("cat.obj");

    Mesh cat(Vector(21,-10,0), 0.6, triangleMesh.vertices, triangleMesh.normals, triangleMesh.uvs, triangleMesh.indices, triangleMesh.vertexcolors);

    Camera camera(Vector(0,0,55), H, W, 1.0472);

    Light light(Vector(-10, 20, 40), 2e10);

    Sphere ball(Vector(21,0,0), 10, Vector(1,1,1), false, false);
    Sphere reflection_ball(Vector(-21,0,0), 10, Vector(1,1,1), true, false);
    Sphere refraction_ball(Vector(0,0,0), 10, Vector(1,0,0), false, true);
    Sphere red(Vector(0,1000,0), 940, Vector(1,0,0), false, false);
    Sphere green(Vector(0,0,-1000), 940, Vector(0,1,0), false, false);
    Sphere yellow(Vector(1000,0,0), 940, Vector(1,1,0), false, false);
    Sphere blue(Vector(0,-1000,0), 990, Vector(0,0,1), false, false);
    Sphere pink(Vector(0,0,1000), 940, Vector(1,0.4,0.7), false, false);
    Sphere lightblue(Vector(-1000,0,0), 940, Vector(0.7,0.85,0.9), false, false);

    std::vector<Geometry *> arrayS;

    if (scene_nbr == 0){
        arrayS.push_back(&ball);
        arrayS.push_back(&reflection_ball);
        arrayS.push_back(&refraction_ball);
    }
    else if (scene_nbr == 1) {arrayS.push_back(&cat);}
    else {arrayS.push_back(&ball);}

    arrayS.push_back(&red);
    arrayS.push_back(&green);
    arrayS.push_back(&yellow);
    arrayS.push_back(&blue);
    arrayS.push_back(&pink);
    arrayS.push_back(&lightblue);

    Scene scene(arrayS, light);

    int max_path_length = 5;
    int n_rays = 10;
    double gamma = 2.2;

    //#pragma omp parallel for schedule(dynamic,1) 
    for (int i=0; i<H; i++) {
        cout << i << endl;
        for (int j=0; j<W; j++) {
            double x,y;
            boxMuller(0.4,x,y);
            Vector color;
            for (int k = 0; k < n_rays; k++) {
                Vector d = Vector(camera.Q[0]+ x + j+ 0.5 - (double) camera.W/2, camera.Q[1] - i - y - 0.5 + (double) camera.H/2, camera.Q[2] - (double) camera.W/(2*tan( (double) camera.alpha/2)));
                Vector u = divide(d-camera.Q, norm(d-camera.Q));
                Ray ray = Ray(camera.Q, u);
                color += scene.getColor(ray , max_path_length);
            }
            image[(i*W + j) * 3 + 0] = min(255., pow(color[0]/n_rays, 1/gamma)); 
            image[(i*W + j) * 3 + 1] = min(255., pow(color[1]/n_rays, 1/gamma)); 
            image[(i*W + j) * 3 + 2] = min(255., pow(color[2]/n_rays, 1/gamma)); 
        }
    }
 
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
 
    return 0;
}
