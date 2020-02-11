//
//  Skeleton3D.hpp
//  projetc
//
//  Created by Coralie DAVID on 27/10/2017.
//  Copyright © 2017 Manon Romain. All rights reserved.
//
#pragma once
#ifndef Skeleton3D_hpp
#define Skeleton3D_hpp

#include <vector>

#include "Quaternion.hpp"
#include "Joint3D.hpp"
#include "Skeleton2D.hpp"

using namespace std;
//typedef CGAL::Simple_cartesian<double>  Kernel;
//typedef Kernel::Point_3 Point_3;
//typedef Kernel::Vector_3 Vector_3;

class Skeleton3D{
    Joint3D root;
    vector<Skeleton3D> children;
    String name;
    public :
    
        Skeleton3D(Joint3D& r, vector<Skeleton3D>& c, String n=""){
            root = r;
            children = c;
            name = n;
        }
        void add_child(Skeleton3D& c){
            this->children.push_back(c);
        }
        size_t get_children_size(){
            return children.size();
        }
        Skeleton3D* get_child(int i){
            return &children[i];
        }
        String get_name(){
            return name;
        }
        Joint3D* get_root() {
            return &root;
        }
        void rotate(Quaternion q);
        void updateAbsolutePosition(vector<Quaternion> q = vector<Quaternion>(), vector<Point3f> t =vector<Point3f>());
        Skeleton2D project(Mat& d);
        void transform(Mat& h);
        void transform_translate(int x, int y, int z);

};


#endif /* Skeleton3D_hpp */
