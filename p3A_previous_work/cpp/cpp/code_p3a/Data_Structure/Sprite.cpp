//
//  Sprite.cpp
//  projetc
//
//  Created by Coralie DAVID on 04/01/2018.
//  Copyright © 2018 Manon Romain. All rights reserved.
//

#include "Sprite.hpp"

void Sprite::incrementPos(){
    if (full_cycle || up){
        pos = (pos + 1)%number;
    }
}
void Sprite::preload_images(String ext){
    //Load only images selected in path
    for (int i = 0; i<number; i++){
        stringstream query;
        query << image_path << path[i] << ext << ".png";
        Mat im = imread(query.str(), -1);
        //cout<<im.size()<<endl;
        try{
            resize(im, im, SIZE_I);
            preloaded_images.push_back(im);
        }
        catch (Exception e){
            cout<<"Couldn't open "<<query.str()<<endl;
            break;
        }
    }
}

void Sprite::display_one_vert(double theta){
    Mat im = preloaded_images[pos];
    Mat M;
    Point2f src_pts[3], dest_pts[3];
    src_pts[0] = Point2f(im.rows/2, im.cols/2);
    dest_pts[0] = Point2f(im.rows/2, im.cols/2);
    src_pts[1] = Point2f(im.rows/2, im.cols/2+5);
    dest_pts[1] = Point2f(im.rows/2, im.cols/2+5);
    src_pts[2] = Point2f(im.rows/2+1, im.cols/2);
    dest_pts[2] = Point2f(im.rows/2+cos(theta*PI/180), im.cols/2+sin(theta*PI/180));
    
    M = getAffineTransform(src_pts, dest_pts);

    Mat out;
    imshow("in", im);
    warpAffine(im, out, M, SIZE_I);
    imshow("out", out);
    incrementPos();
}

void Sprite::display_one_horiz(double theta){
    Mat im = preloaded_images[pos];
    Mat M;
    Point2f src_pts[3], dest_pts[3];
    src_pts[0] = Point2f(im.rows/2, im.cols/2);
    dest_pts[0] = Point2f(im.rows/2, im.cols/2);
    src_pts[1] = Point2f(im.rows/2+5, im.cols/2);
    dest_pts[1] = Point2f(im.rows/2+5, im.cols/2);
    src_pts[2] = Point2f(im.rows/2, im.cols/2+1);
    dest_pts[2] = Point2f(im.rows/2+cos(theta*PI/180), im.cols/2+sin(theta*PI/180));
    
    M = getAffineTransform(src_pts, dest_pts);
    
    Mat out;
    warpAffine(im, out, M, SIZE_I);
    imshow("horiz", out);
}
