#include <iostream>
#include <vector>
#include <math.h> 

#include "reduce_image.cpp"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

using namespace std;

int main() {

    string file = "pitieJones.png";

    int W, H, c;

    unsigned char *data = stbi_load(file.c_str(), &W, &H, &c, 3); 

    vector<unsigned char> image(W*H*3,0);
    int i = 0;
    for (int w=0; w<W; w++){
        for (int h=0; h<H; h++){
            image[(w*W + h) * 3 + 0] = data[i];
            image[(w*W + h) * 3 + 1] = data[i+1];
            image[(w*W + h) * 3 + 2] = data[i+2];
            i += 3; 
        }
    }

    int nbr_crop = 10;

    vector<unsigned char> image_new = reduce_image(image, H, W, nbr_crop);

    stbi_write_png("image.png", W, H, 3, &image_new[0], 0);

    return 0;
}
