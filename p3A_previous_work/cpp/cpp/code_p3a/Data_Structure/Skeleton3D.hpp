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
#include <map>

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
    vector<int> hierarchy;      // give the position w.r.t the root of the skeleton
    float min_angle;       // minimal angle in degree with the horizontal line
    float max_angle;       // \in [-180,180]
    public :
    
        Skeleton3D(Joint3D& r, vector<Skeleton3D>& c, String n){
            root = r;
            children = c;
            name = n;
            min_angle = 0;
            max_angle = 0;
        }
        void add_child(Skeleton3D& c){
            this->children.push_back(c);
        }
        void add_node(int d) {
            this->hierarchy.push_back(d);
        }
        void add_constraint(float min, float max) {
            this->min_angle = min;
            this->max_angle = max;
        }
        vector<float> get_angle_constraints() {
            vector<float> v(min_angle, max_angle);
            return v;
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
        size_t get_hierarchy_size() {
            return hierarchy.size();
        }
        int get_node(int i) {
            return hierarchy[i];
        }
        vector<int> get_hierarchy() {
            return hierarchy;
        }
        void rotate(Quaternion q);
        void updateAbsolutePosition(vector<Quaternion> q = vector<Quaternion>(), vector<Point3f> t =vector<Point3f>());
        Skeleton2D project();
        Skeleton2D project_individual();
        void transform(Mat& h);
        void transform_translate(int x, int y, int z);
        void create_hierarchy();
        void add_angle_constraints();
};


#endif /* Skeleton3D_hpp */
