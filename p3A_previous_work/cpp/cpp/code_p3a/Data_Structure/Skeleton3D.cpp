//
//  Skeleton3D.cpp
//  projetc
//
//  Created by Coralie DAVID on 27/10/2017.
//  Copyright © 2017 Manon Romain. All rights reserved.
//

#include "Skeleton3D.hpp"
#include <map>
#include "../constantes.hpp"

void Skeleton3D::transform(Mat& h){
    Point3f p = h*(root.getPos());
    if (name != "hip") {
        root.setPos(p);
    }
    for (size_t i=0; i<children.size(); i++){
        children[i].transform(h);
    }
    updateAbsolutePosition();
}

void Skeleton3D::transform_translate(int x, int y, int z) {
    Point3f h = Point3f(x,y,z);
    Point3f p = h + (root.getPos());
    if (name == "hip") {        // only need to move the hip, as everything is attached to it
        root.setPos(p);
    }
    for (size_t i = 0; i < children.size(); i++) {
        children[i].transform_translate(x,y,z);
    }
    updateAbsolutePosition();
}

void Skeleton3D::create_hierarchy() {
    // also add angles constraints
    for (size_t i = 0; i < children.size(); i++) {
        children[i].hierarchy = hierarchy;  // hierarchy of the parent
        children[i].add_node(i);            // add its own hierarchy w.r.t the parent
        children[i].create_hierarchy();     // iterate
    }
}

void Skeleton3D::add_angle_constraints(){
    if (min_angle == max_angle) {                       // if no constraint on angles
        float alpha = get_angle(root.getPos()) * rtd;    // get angle in degrees
        add_constraint(alpha - 45, alpha + 45);         // +- 35 degrees wrt initial position
        //cout << name << ": " << min_angle << "; " << alpha << "; " << max_angle << endl;
    }
    //angles = get_angle_constraints();
    for (size_t i = 0; i < children.size(); i++) {
        children[i].add_angle_constraints();
    }
}



void Skeleton3D::updateAbsolutePosition(vector<Quaternion> q, vector<Point3f> t){
    Point3f rel_pos = this->root.getPos();
    Point3f absolute_position = rel_pos;
    Quaternion curr_q =this->root.getQuat();

    Point3f r;

    //cout << "t: " << t << endl;
    for (int i = q.size()-1; i>=0; i--){
        //Quaternion::applyQuaternion(q[i], absolute_position);       // f(abs) = q[i] * abs * q[i]^-1  => rotation around q[i] axis
        absolute_position += t[i];
    }
    this->root.setAbsPos(absolute_position);
    q.push_back(curr_q);
    t.push_back(rel_pos);
    
    for (size_t i=0; i<children.size(); i++){
        children[i].updateAbsolutePosition(q, t);
    }
};

Skeleton2D Skeleton3D::project(){
    vector<Skeleton2D> liste;
    for (size_t i=0; i<children.size(); i++){
        liste.push_back(children[i].project());
    }
    Point3f aux = this->root.getAbsPos();
    Point3f p_project(aux.x, aux.y, 1);
    return Skeleton2D(p_project, liste, name);
};

Skeleton2D Skeleton3D::project_individual() {
    // only project 1 bone, not the entire skeleton
    vector<Skeleton2D> liste;
    vector<Skeleton2D> liste_temp;
    Point3f aux_temp = this->root.getAbsPos();
    Point3f project_temp(aux_temp.x, aux_temp.y, 1);
    liste.push_back(Skeleton2D(project_temp, liste_temp, name));
    Point3f aux = this->root.getAbsPos() - this->root.getPos();     // <=> get Abs Position of the parent
    Point3f p_project(aux.x, aux.y, 1);
    return Skeleton2D(p_project, liste, name);
}

void Skeleton3D::rotate(Quaternion q){
    Quaternion qf;
    cout << "root and q: " << this->name << "; " << q << endl;
    Quaternion::multiplyQuaternion(this->root.getQuat(), q, qf);
    cout << "Debug rotate: " << this->root.getQuat() << "; " << qf;
    this->root.setQuat(qf);
    cout << "; " << this->root.getQuat() << endl;
};

