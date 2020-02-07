//
//  Skeleton3D.cpp
//  projetc
//
//  Created by Coralie DAVID on 27/10/2017.
//  Copyright © 2017 Manon Romain. All rights reserved.
//

#include "Skeleton3D.hpp"

void Skeleton3D::transform(Mat& h){
    //cout << "DEBUG h, root: " << name << "\n" << h << "\n" << root.getPos() << endl;
    Point3f p = h*(root.getPos());
    //cout << "multiplication (h*root): " << p << endl;
    if (name != "hip") {
        root.setPos(p);
    }
    for (size_t i=0; i<children.size(); i++){
        children[i].transform(h);
    }
    //cout << endl << name << endl;
    updateAbsolutePosition();
}

void Skeleton3D::transform_translate(int x, int y, int z) {
    Point3f h = Point3f(x,y,z);
    //cout << "DEBUG h: " << name << root.getPos() << endl;
    Point3f p = h + (root.getPos());
    //Point3f p = h;
    if (name == "hip") {        // only need to move the hip, as everything is attached to it
        root.setPos(p);
    }
    //cout << "root pos: " << root.getPos() << endl;
    for (size_t i = 0; i < children.size(); i++) {
        children[i].transform_translate(x,y,z);
    }
    updateAbsolutePosition();
}

void Skeleton3D::updateAbsolutePosition(vector<Quaternion> q, vector<Point3f> t){
    Point3f rel_pos = this->root.getPos();
    Point3f absolute_position = rel_pos;
    Quaternion curr_q =this->root.getQuat();
    //cout << "update Absolute position: " << q.size() << t << endl;

    Point3f r;

    //cout << "t: " << t << endl;
    for (int i = q.size()-1; i>=0; i--){
        //Quaternion::applyQuaternion(q[i], absolute_position);       // f(abs) = q[i] * abs * q[i]^-1  => rotation around q[i] axis
        absolute_position += t[i];
        //cout << name << " absolute position: " << absolute_position << endl;
    }
    this->root.setAbsPos(absolute_position);
    q.push_back(curr_q);
    t.push_back(rel_pos);
    
    for (size_t i=0; i<children.size(); i++){
        children[i].updateAbsolutePosition(q, t);
    }
};

Skeleton2D Skeleton3D::project(Mat& h){
    vector<Skeleton2D> liste;
    for (size_t i=0; i<children.size(); i++){
        liste.push_back(children[i].project(h));
    }
    //cout << "abs pos 2D skeleton: " << this->name << this->root.getAbsPos() << endl;
    Point3f aux = h*this->root.getAbsPos();
    Point3f p_project(aux.x, aux.y, 1);
    return Skeleton2D(p_project, liste, name);
};

void Skeleton3D::rotate(Quaternion q){
    Quaternion qf;
    Quaternion::multiplyQuaternion(this->root.getQuat(), q, qf);
    this->root.setQuat(qf);
    
};

