#include <iostream>
#include <vector>
#include <math.h> 

#include "clip.cpp"
#include "voronoi.cpp"

using namespace std;

int main(){
        vector<Vector> pointSubj;
        pointSubj.push_back(Vector(0.0, 0.1, 0));
        pointSubj.push_back(Vector(0.1, 0.2, 0));
        pointSubj.push_back(Vector(0.3, 0, 0));
        pointSubj.push_back(Vector(0.5, 0.1, 0));
        pointSubj.push_back(Vector(0.4, 0.4, 0));
        pointSubj.push_back(Vector(0.1, 0.3, 0));

        Polygon subjectPolygon = Polygon(pointSubj);

        vector<Vector> pointClip;
        pointClip.push_back(Vector(0.1, 0.1, 0));
        pointClip.push_back(Vector(0.1, 0.4, 0));
        pointClip.push_back(Vector(0.4, 0.4, 0));
        pointClip.push_back(Vector(0.4, 0.1, 0));

        Polygon clipPolygon = Polygon(pointClip);

        vector<Polygon> polygons;
        polygons.push_back(subjectPolygon);
        polygons.push_back(clipPolygon);

        vector<Vector> empty;
        save_svg(polygons, "before.svg", empty);

        vector<Polygon> poly;
        Polygon newPoly = clip(subjectPolygon, clipPolygon);
        poly.push_back(newPoly);

        save_svg(poly, "after.svg", empty);

        uint n = 50;
        vector<Vector> pointBox;
        pointBox.push_back(Vector(0.1, 0.1));
        pointBox.push_back(Vector(0.1, 0.8));
        pointBox.push_back(Vector(0.8, 0.8));
        pointBox.push_back(Vector(0.8, 0.1));
        Polygon box(pointBox);

        
        vector<Vector> p(n);
        for (int i = 0; i < n; i++) {
            double x = (double) rand() * 0.7 / RAND_MAX + 0.1;
            double y = (double) rand() * 0.7 / RAND_MAX + 0.1;
            p[i] = Vector(x, y);
        }
        Polygon points(p);
        save_svg(voronoi_diagram(points, box), "diagram.svg", p);

        vector<double> weights(n);
        for (int i = 0; i < n; i++) {
            if (p[i][0] < 0.2 || p[i][0] > 0.7 || p[i][1] < 0.2 || p[i][1] > 0.7) {
                weights[i] = 0.8;
            } else {
                weights[i] = 1;
            }
        }
        save_svg(voronoi_diagram(points, box, weights), "power_diagram.svg", p);

    return 0;
}
