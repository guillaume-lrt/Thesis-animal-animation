//
//  Quaternion.cpp
//  projetc
//
//  Created by Coralie DAVID on 27/10/2017.
//  Copyright © 2017 Manon Romain. All rights reserved.
//

#include "Quaternion.hpp"
#include <iostream>

double Quaternion::norm2(){
    double x = this->x[0];
    double y = this->x[1];
    double z = this->x[2];
    double r = this->x[3];
    return x*x+y*y+z*z+r*r;
}
Quaternion Quaternion::inv(){
    double n2 = this->norm2();
    return Quaternion((this->x[0])/n2, -(this->x[1])/n2, -(this->x[2])/n2, -(this->x[3])/n2);
}
void Quaternion::multiplyQuaternion(const Quaternion& q1,const Quaternion& q2, Quaternion& q)
{
    // First quaternion q1 (x1 y1 z1 r1)
    const float r1 = q1(0);
    const float x1 = q1(1);
    const float y1 = q1(2);
    const float z1 = q1(3);
    
    // Second quaternion q2 (x2 y2 z2 r2)
    const float r2 = q2(0);
    const float x2 = q2(1);
    const float y2 = q2(2);
    const float z2 = q2(3);
    
    q(0) = r1*r2 - x1*x2 - y1*y2 - z1*z2;   // r component
    q(1) = x1*r2 + r1*x2 + y1*z2 - z1*y2;   // x component
    q(2) = r1*y2 - x1*z2 + y1*r2 + z1*x2;   // y component
    q(3) = r1*z2 + x1*y2 - y1*x2 + z1*r2;   // z component
    
}
void Quaternion::applyQuaternion(Quaternion& q, Point3f& p){
    // First quaternion q1 (x1 y1 z1 r1)
    cout << "quaternion: " << q << " and point (absolute position): " << p << endl;
    Quaternion aux(p);      // convert p into a quaternion
    Quaternion a1, a2;      // initialize as (0,0,0,0)
    Quaternion::multiplyQuaternion(q, aux, a1);         // a1 = q * aux
    Quaternion::multiplyQuaternion(a1, q.inv(), a2);    // a2 = a1 * q^-1
    p = Point3f(a2(1), a2(2), a2(3));   // p = q * quad(p) * q^-1

    cout << "new point: " << p << endl;
    return ;
}
