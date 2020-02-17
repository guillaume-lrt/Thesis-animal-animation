//
//  Skeleton2D.hpp
//  projetc
//
//  Created by Coralie DAVID on 27/10/2017.
//  Copyright © 2017 Manon Romain. All rights reserved.
//
//#pragma once
#ifndef Skeleton2D_hpp
#define Skeleton2D_hpp

#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>

#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/videoio/videoio.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "utils.hpp"

using namespace std;
using namespace cv;
class Skeleton2D{
    Point3f pos;
    vector<Skeleton2D> children;
    String name;
    vector<double> minMax;
    void auxMat(int width, int height, Mat& dest, Mat& d, string shape);
    void oneMat(int width, int height, Mat& dest, Mat& debug_im, int i, string shape);
    public :
    Skeleton2D(Point3f& p, vector<Skeleton2D> c, String n){
        this->pos = p/p.z;
        this->children = c;
        this->name = n;
        double minx, maxx, miny, maxy;
        minx = p.x; maxx = p.x;
        miny = p.y; maxy = p.y;
        for (size_t i=0; i<children.size(); i++){
            minx = min(minx, children[i].minMax[0]);
            maxx = max(maxx, children[i].minMax[1]);
            miny = min(miny, children[i].minMax[2]);
            maxy = max(maxy, children[i].minMax[3]);
        }
        this->minMax.push_back(minx);
        this->minMax.push_back(maxx);
        this->minMax.push_back(miny);
        this->minMax.push_back(maxy);
    }
    void add_child(Skeleton2D& s){
        this->children.push_back(s);
        this->minMax[0] = min(this->minMax[0], s.minMax[0]);
        this->minMax[1] = max(this->minMax[1], s.minMax[1]);
        this->minMax[2] = min(this->minMax[2], s.minMax[2]);
        this->minMax[3] = max(this->minMax[3], s.minMax[3]);
        
    }
    void setPos(Point3f& p){
        this->pos = p;
    }
    String get_name(){
        return name;
    }
    void updateMinMax();
    void transform(Mat& h);
    void normalize(int width, int height);
    Mat toMat(int width, int height, bool show = false, string shape = "ellipse");
    
    
};
#endif /* Skeleton2D_hpp */
