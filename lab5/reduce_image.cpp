#include <iostream>
#include <vector>
#include <math.h> 

#include "vector.cpp"

using namespace std;

double I(int x, int y, int& W, int& H, vector<unsigned char>& image){
    double R = image[(x*W + y) * 3 + 0];
    double G = image[(x*W + y) * 3 + 1];
    double B = image[(x*W + y) * 3 + 2];
    return R+G+B;
}

vector<double> Energy_map(int& W, int& H, vector<unsigned char>& image){
    vector<double> energy_map(W*H, 0);

    for (int i=0; i<H; i++){
        for (int j=0; j<W; j++){
            if (i >= 1 && i < H-1) energy_map[i*W+j] += abs(I(i-1, j, W, H, image) - I(i+1, j, W, H, image));
            else if (i == 0) energy_map[i*W+j] += abs(I(i, j, W, H, image) - I(i+1, j, W, H, image));
            else if (i == H-1) energy_map[i*W+j] += abs(I(i, j, W, H, image) - I(i-1, j, W, H, image));

            if (j >= 1 && j < W-1) energy_map[i*W+j] += abs(I(i, j-1, W, H, image) - I(i, j+1, W, H, image));
            else if (j == 0) energy_map[i*W+j] += abs(I(i, j, W, H, image) - I(i, j+1, W, H, image));
            else if (j == W-1) energy_map[i*W+j] += abs(I(i, j, W, H, image) - I(i, j-1, W, H, image));
        }
    }

    return energy_map;
}

vector<double> cumulated(int& W, int& H, vector<unsigned char>& image){
    vector<double> energy_map = Energy_map(W, H, image);
    vector<double> cumulated_energy_map(W*H, 0);

    for (int i=0; i<H; i++){
        for (int j=0; j<W; j++){
            if (i==0) cumulated_energy_map[j] = energy_map[j];
            else{
                if (j==0) cumulated_energy_map[i*W+j] = min(cumulated_energy_map[(i-1)*W+j], cumulated_energy_map[(i-1)*W+j+1]) + energy_map[i*W+j];
                else if (j==W-1) cumulated_energy_map[i*W+j] = min(cumulated_energy_map[(i-1)*W+j-1], cumulated_energy_map[(i-1)*W+j]) + energy_map[i*W+j];
                else cumulated_energy_map[i*W+j] = min(cumulated_energy_map[(i-1)*W+j-1], min(cumulated_energy_map[(i-1)*W+j], cumulated_energy_map[(i-1)*W+j+1])) + energy_map[i*W+j];
            }    
        }
    }
    return cumulated_energy_map;
}

vector<int> find_seam(int& W, int& H, vector<unsigned char>& image){
    vector<double> cumulated_energy_map = cumulated(W, H, image);
    double min = numeric_limits<double >::max();
    int min_index = 0;
    vector<int> seam;

    for (int j=0; j<H; j++){
        if (cumulated_energy_map[(W-1)*W + j] < min){
            min = cumulated_energy_map[j];
            min_index = j;
        }
    }

    seam.insert(seam.begin(), min_index);

    for (int i=H-1; i>0; i--){
        int y = seam[0];

        if (y==0){ 
            if (cumulated_energy_map[(i-1)*W + y] < cumulated_energy_map[(i-1)*W + y+1]){
                seam.insert(seam.begin(), y);
            } else {seam.insert(seam.begin(), y+1); }
        }
        else if (y==W-1){ 
            if(cumulated_energy_map[(i-1)*W + y-1] < cumulated_energy_map[(i-1)*W + y]){
                seam.insert(seam.begin(), y-1);
            } else {seam.insert(seam.begin(), y); }
        }
        else { 
            double val = numeric_limits<double >::max();
            int index = 0;
            for(int k=0; k<3; k++){
                if (cumulated_energy_map[(i-1)*W + y-1+k] < val){
                    val = cumulated_energy_map[(i-1)*W + y-1+k];
                    index = k;
                }
            }
            seam.insert(seam.begin(), y-1+index);
        }

    }
    return seam;
}

vector<unsigned char> reduce_image(vector<unsigned char>& image, int& H, int& W, int nbr){

    cout << nbr << endl;

    vector<unsigned char> new_image((W-nbr)*H*3,0);

    for(int k=0; k<nbr; k++){
        vector<int> seam = find_seam(W, H, image);

        for (int i=0; i<H; i++){
            int remove = seam[i];
            for (int j=0; j<W; j++){
                if (j < remove){
                    new_image[(i*(W-1)+j)*3 + 0] = image[(i*(W)+j)*3 + 0]; 
                    new_image[(i*(W-1)+j)*3 + 1] = image[(i*(W)+j)*3 + 1];
                    new_image[(i*(W-1)+j)*3 + 2] = image[(i*(W)+j)*3 + 2];  
                }
                else{
                    new_image[(i*(W-1)+j)*3 + 0] = image[(i*(W)+j+1)*3 + 0]; 
                    new_image[(i*(W-1)+j)*3 + 1] = image[(i*(W)+j+1)*3 + 1];
                    new_image[(i*(W-1)+j)*3 + 2] = image[(i*(W)+j+1)*3 + 2];  
                }
            }
        }
    }
    return new_image;
}
