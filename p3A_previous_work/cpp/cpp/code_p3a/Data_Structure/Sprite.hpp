//
//  Sprite.hpp
//  projetc
//
//  Created by Coralie DAVID on 04/01/2018.
//  Copyright © 2018 Manon Romain. All rights reserved.
//

#ifndef Sprite_hpp
#define Sprite_hpp

#include <stdio.h>
#include <vector>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/videoio/videoio.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "../constantes.hpp"
#include "utils.hpp"

using namespace std;
using namespace cv;
class Sprite{
    Vec3f v_n;
    String image_path;
    vector<int> path;
    int pos;
    int number;
    bool full_cycle, up;
    //Preloaded images
    vector<Mat> preloaded_images;
    

    void incrementPos();
    void preload_images(String ext = "");
    
    public:
        void display_one_horiz(double theta);
        void display_one_vert(double theta);
        Sprite(String n, vector<int> p, bool cycle, Vec3f vec_n){
            image_path = n;
            path = p;
            number = p.size();
            pos = 0; up = true;
            full_cycle = cycle;
            
            v_n = vec_n;
            
            preloaded_images = vector<Mat>();
            preload_images();
            cout<<"Preloaded images successfully"<<endl;
        }
    
    
};
#endif /* Sprite_hpp */
