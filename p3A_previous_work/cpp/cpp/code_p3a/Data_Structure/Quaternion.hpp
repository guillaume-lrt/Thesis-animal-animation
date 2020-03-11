//
//  Quaternion.hpp
//  projetc
//
//  Created by Coralie DAVID on 27/10/2017.
//  Copyright © 2017 Manon Romain. All rights reserved.
//
//#pragma once
#ifndef Quaternion_hpp
#define Quaternion_hpp

#include <vector>
//#include <CGAL/Simple_cartesian.h>
//
//typedef CGAL::Simple_cartesian<double>  Kernel;
//typedef Kernel::Point_3 Point_3;

#include <opencv2/imgcodecs.hpp>
#include <opencv2/videoio/videoio.hpp>
#include <opencv2/highgui/highgui.hpp>


using namespace std;
using namespace cv;

class Quaternion{
    vector<float> x;
    public:
        Quaternion(){
            this->x.push_back(0);
            this->x.push_back(0);
            this->x.push_back(0);
            this->x.push_back(0);
        };
        Quaternion(float xa, float xb, float xc, float xd){
            this->x.push_back(xa);
            this->x.push_back(xb);
            this->x.push_back(xc);
            this->x.push_back(xd);
        };
        Quaternion(Point3f& p){
            this->x.push_back(0);
            this->x.push_back(p.x);
            this->x.push_back(p.y);
            this->x.push_back(p.z);
        };
        Quaternion(double angle, Vec3f& v){
            this->x.push_back(float(cos(angle)));
            this->x.push_back(float(sin(angle)*v[0]));
            this->x.push_back(float(sin(angle)*v[1]));
            this->x.push_back(float(sin(angle)*v[2]));
        };
        std::vector<float> getX(){
            return this->x;
        }
        double norm2();
        Quaternion inv();
        inline float operator()(int i) const { return x[i]; }
        inline float& operator()(int i) { return x[i]; }
        friend ostream& operator<<(ostream& os, const Quaternion& q){
            os<<"q="<<q.x[0]<<", "<<q.x[1]<<", "<<q.x[2]<<", "<<q.x[3];
            return os;
        }
    
        static void applyQuaternion(Quaternion& q, Point3f& p);
        static void multiplyQuaternion(const Quaternion& q1,const Quaternion& q2, Quaternion& q);

};
#endif /* Quaternion_hpp */
