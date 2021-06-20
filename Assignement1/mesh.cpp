#pragma once

#include <iostream>  
#include <vector>
#include <cmath>
#include <math.h>
#include <list>

#include "mesh_loader.cpp"
#include "boundingbox.cpp"

using namespace std;

class Node {
    public:
        Node *child_left;
        Node *child_right;
        BoundingBox bbox;
        int starting_triangle;
        int ending_triangle;

        Node(){
            this->child_left = NULL;
            this->child_right = NULL;
        }
};

class Mesh : public Geometry{
    public:

        Mesh(Vector T, double scale, vector<Vector> vertices, vector<Vector> normals, vector<Vector> uvs, vector<TriangleIndices> triangles, std::vector<Vector> vertexcolors){
            this->vertices = vertices;
            this->normals = normals;
            this->uvs = uvs;
            this->triangles = triangles;
            this->albedo = Vector(4.,4.,4.);
            this->T = T;
            this->scale = scale;
            int start = 0;
            int end = triangles.size();
            root = new Node();
            box = compute_bbox(start, end, T, scale);
            build_bvh(root, start, end);
        }

        void build_bvh(Node *node, int& starting_triangle, int& ending_triangle) {
            node->bbox = compute_bbox(starting_triangle, ending_triangle,T, scale);
            node->starting_triangle = starting_triangle;
            node->ending_triangle = ending_triangle;
            Vector diag = compute_diag(node->bbox);
            Vector middle_diag = node->bbox.min + diag * 0.5;
            int longest_axis = get_longest(diag);
            int pivot_index = starting_triangle;
            for (int i = starting_triangle; i < ending_triangle; i++) {
                Vector barycenter = compute_barycenter(i, scale, T);
                if (barycenter[longest_axis] < middle_diag[longest_axis]) {
                    swap(triangles[i], triangles[pivot_index]);
                    pivot_index++;
                }
            }
            if ( pivot_index <= starting_triangle || pivot_index >= ending_triangle - 1 || ending_triangle - starting_triangle < 5){return;}
            node->child_left = new Node();
            node->child_right = new Node();
            build_bvh(node->child_left, starting_triangle, pivot_index);
            build_bvh(node->child_right, pivot_index, ending_triangle);
        }
    
        Intersection intersect_old(Ray ray){
            Intersection inter;
            inter.result = false;
            double min_distance = numeric_limits<double >::max();
            double temp;
            if (this->box.intersect(ray, temp)){
                for(int i=0; i<triangles.size(); i++){
                    Vector A = vertices[triangles[i].vtxi] * scale + T;
                    Vector B = vertices[triangles[i].vtxj] * scale + T;
                    Vector C = vertices[triangles[i].vtxk] * scale + T;

                    Vector e1 = B-A;
                    Vector e2 = C-A;
                    Vector N = cross(e1,e2);

                    Vector AOu = cross(A-ray.O, ray.u);
                    double uN = dot(ray.u,N);

                    double beta = dot(e2, AOu)/uN;
                    double gamma = - dot(e1, AOu)/uN;
                    double alpha = 1. - beta - gamma;

                    if (alpha >= 0. && beta >= 0. && gamma >= 0.) {
                        double t = dot(A-ray.O, N)/uN;

                        if (t < min_distance && t > 0.){ 
                            min_distance = t;
                            inter = Intersection(true, ray.O + ray.u*t, divide(N,norm(N)), t, albedo);
                        }
                    }
                }
            }
            return inter; 
        }

        Intersection intersect(Ray ray){
            Intersection inter(false);
            double temp;
            if (!root->bbox.intersect(ray, temp)){return Intersection(false);}
            list<Node*> nodes_to_visit;
            nodes_to_visit.push_front(root);
            Node* curNode;
            double best_inter_distance = std::numeric_limits<double >::max();
            while(!nodes_to_visit.empty()){
                curNode = nodes_to_visit.back();
                nodes_to_visit.pop_back();
                double inter_distance;
                if (curNode->child_left){
                    if(curNode->child_left->bbox.intersect(ray, inter_distance)){
                        if (inter_distance < best_inter_distance ){
                            nodes_to_visit.push_back(curNode->child_left);
                        }
                    }
                    if(curNode->child_right->bbox.intersect(ray, inter_distance)){
                        if (inter_distance < best_inter_distance ){
                            nodes_to_visit.push_back(curNode->child_right);
                        }
                    }
                }
                else{
                    double min_distance = std::numeric_limits<double>::max();
                    for(int i=curNode->starting_triangle; i<curNode->ending_triangle; i++){
                        Vector A = vertices[triangles[i].vtxi] * scale + T;
                        Vector B = vertices[triangles[i].vtxj] * scale + T;
                        Vector C = vertices[triangles[i].vtxk] * scale + T;

                        Vector e1 = B-A;
                        Vector e2 = C-A;
                        Vector N = cross(e1,e2);

                        Vector AOu = cross(A-ray.O, ray.u);
                        double uN = dot(ray.u,N);

                        double beta = dot(e2, AOu)/uN;
                        double gamma = - dot(e1, AOu)/uN;
                        double alpha = 1. - beta - gamma;

                        if (alpha >= 0. && beta >= 0. && gamma >= 0.) {
                            double t = dot(A-ray.O, N)/uN;

                            if (t < min_distance && t > 0.){ 
                                min_distance = t;
                                inter = Intersection(true, ray.O + ray.u*t, divide(N,norm(N)), t, albedo);
                            }
                        }
                    }
                    if (min_distance < best_inter_distance){
                        best_inter_distance = min_distance;
                        inter.t = best_inter_distance;
                    }
                }
            }
            return inter;
        }

    private:
        double scale;
        Vector T;
        std::vector<Vector> vertices;
        std::vector<Vector> normals;
        std::vector<Vector> uvs;
        std::vector<TriangleIndices> triangles;
        BoundingBox box;
        Node *root;

        int get_longest(Vector &diag){
            int longest = 0;
            double max = abs(diag[0]);
            for (int i = 1; i < 3; i++) {
                if (abs(diag[i]) > max) {
                max = abs(diag[i]);
                longest = i;
                }
            }
            return longest;
        }

        BoundingBox compute_bbox(int& starting_triangle, int& ending_triangle, Vector T, double scale){
            BoundingBox bbox;
            Vector max = Vector(-numeric_limits<double >::max(), -numeric_limits<double >::max(), -numeric_limits<double >::max());
            Vector min = Vector(numeric_limits<double >::max(), numeric_limits<double >::max(), numeric_limits<double >::max());
            for(int i=starting_triangle; i<ending_triangle; i++){
                vector<Vector> tempList;
                tempList.push_back(vertices[triangles[i].vtxi]);
                tempList.push_back(vertices[triangles[i].vtxj]);
                tempList.push_back(vertices[triangles[i].vtxk]);
                for (int j=0; j<3; j++){
                    if (tempList[j][0] * scale + T[0] < min[0]) min[0] = tempList[j][0] * scale + T[0];
                    else if (tempList[j][0] * scale + T[0] > max[0]) max[0] = tempList[j][0] * scale + T[0];
                    if (tempList[j][1] * scale + T[1] < min[1]) min[1] = tempList[j][1] * scale + T[1];
                    else if (tempList[j][1] * scale + T[1] > max[1]) max[1] = tempList[j][1] * scale + T[1];
                    if (tempList[j][2] * scale + T[2] < min[2]) min[2] = tempList[j][2] * scale + T[2];
                    else if (tempList[j][2] * scale + T[2] > max[2]) max[2] = tempList[j][2] * scale + T[2];
                }
                
            }
            bbox.min = min;
            bbox.max = max;
            return bbox;
        }

        Vector compute_diag(BoundingBox bbox){
            return bbox.max - bbox.min;
        }

        Vector compute_barycenter(int &i, double &scale, Vector& T){
            Vector vertex1 = vertices[triangles[i].vtxi] * scale + T;
            Vector vertex2 = vertices[triangles[i].vtxj] * scale + T;
            Vector vertex3 = vertices[triangles[i].vtxk] * scale + T;
            Vector barycenter = (vertex1 + vertex2 + vertex3)*(1/3);
            return barycenter;
        }
};
