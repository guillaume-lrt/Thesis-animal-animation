//
//  Joint3D.hpp
//  projetc
//
//  Created by Coralie DAVID on 27/10/2017.
//  Copyright © 2017 Manon Romain. All rights reserved.
//
#pragma once
#ifndef Joint3D_hpp
#define Joint3D_hpp

#include <stdio.h>
#include <iostream>
#include "Quaternion.hpp"
#include <opencv2/imgcodecs.hpp>
#include <opencv2/videoio/videoio.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;
using namespace cv;

class Joint3D{
    Point3f pos;
    Quaternion q;
    Point3f absolute_p;
    public:
        Joint3D(){
            pos = Point3f();
            q = Quaternion();
        };
        Joint3D(Point3f r, Quaternion q);
        Point3f getPos(){
            return this->pos;
        }
        Quaternion getQuat(){
            return this->q;
        }
        Point3f getAbsPos(){
            return this->absolute_p;
        }
        void setAbsPos(Point3f& p){
            this->absolute_p = p;
        }
        void setPos(Point3f& p){
            this->pos = p;
        }
        void setQuat(Quaternion& qf){
            this->q = qf;
        }
};

#endif /* Joint3D_hpp */
