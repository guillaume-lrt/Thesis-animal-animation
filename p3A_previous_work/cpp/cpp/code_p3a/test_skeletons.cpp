//
//  test_skeletons.cpp
//  projetc
//
//  Created by Coralie DAVID on 25/11/2017.
//  Copyright © 2017 Manon Romain. All rights reserved.
//

#include "test_skeletons.hpp"

Skeleton2D test2D(){
    vector<Skeleton2D> vide;
    Point3f ptr(8*10, 5*10, 1);
    Skeleton2D trump(ptr, vide, "trump");
    Point3f n(7.8*10, 2, 1);
    cout << "test skeleton 2D" << endl;
    Skeleton2D nose(n, vide, "nose");
    nose.add_child(trump);
    Point3f h(7.2, 1.5, 1);
    Skeleton2D head(h, vide, "head");
    head.add_child(nose);
    
    Point3f xlf(7, 6, 1);
    Skeleton2D xleft(xlf, vide, "xlf");
    Point3f xrf(5, 6, 1);
    Skeleton2D xright(xrf, vide, "xrf");
    Point3f lf(6.5, 3.6, 1);
    Skeleton2D left(lf, vide, "lf");
    left.add_child(xleft);
    Point3f rf(5.5, 4, 1);
    Skeleton2D right(rf, vide, "rf");
    right.add_child(xright);
    Point3f ne(6, 1.8, 1);
    Skeleton2D neck(ne, vide, "neck");
    neck.add_child(left); neck.add_child(right); neck.add_child(head);
    
    Point3f xlbf(2, 6, 1);
    Skeleton2D xleftback(xlbf, vide, "xlb");
    Point3f xrbf(0, 6, 1);
    Skeleton2D xrightback(xrbf, vide, "xrb");
    Point3f xt(0.5, 3, 1);
    Skeleton2D xtail(xt, vide, "xtail");
    
    Point3f lbf(2.3, 4, 1);
    Skeleton2D leftback(lbf, vide, "lb");
    leftback.add_child(xleftback);
    Point3f rbf(1, 4, 1);
    Skeleton2D rightback(rbf, vide, "rb");
    rightback.add_child(xrightback);
    Point3f t(0.8, 1.8, 1);
    Skeleton2D tail(t, vide, "tail");
    tail.add_child(xtail);
    
    Point3f hip_p(2, 2, 1);
    Skeleton2D hip(hip_p, vide, "hip");
    hip.add_child(neck);hip.add_child(leftback);hip.add_child(rightback);hip.add_child(tail);
    return hip;
}

Skeleton3D test3D(){
    vector<Skeleton3D> vide;
    Quaternion un(1, 0, 0, 0);

    // head part
    Point3f pt_ptr(0, 2, 0);
    Joint3D ptr(pt_ptr, un);
    Skeleton3D trump(ptr, vide,"trump");
    Point3f pt_n(1, 2, 0);
    Joint3D n(pt_n, un);
    Skeleton3D nose(n, vide, "nose");
    nose.add_child(trump);
    Point3f pt_h(2.2, 0.5, 0);
    Joint3D h(pt_h, un);
    Skeleton3D head(h, vide, "head");
    head.add_child(nose);
    
    // front legs
    Point3f pt_xxlf(0.1, 1, 0);
    Joint3D xxlf(pt_xxlf, un);
    Skeleton3D xxleftfront(xxlf, vide, "xxlf");
    Point3f pt_xxrf(-0.1, 1, 0);
    Joint3D xxrf(pt_xxrf, un);
    Skeleton3D xxrightfront(xxrf, vide, "xxrf");
    xxrightfront.add_constraint(45,170);
    Point3f pt_xlf(0.5, 1.4, 0);
    Joint3D xlf(pt_xlf, un);
    Skeleton3D xleftfront(xlf, vide, "xlf");
    xleftfront.add_child(xxleftfront);
    Point3f pt_xrf(-0.5, 1.4, 0);
    Joint3D xrf(pt_xrf, un);
    Skeleton3D xrightfront(xrf, vide, "xrf");
    xrightfront.add_child(xxrightfront);
    Point3f pt_lf(.5, 2.7, -0.8);
    Joint3D lf(pt_lf, un);
    Skeleton3D leftfront(lf, vide, "lf");
    leftfront.add_child(xleftfront);
    Point3f pt_rf(-.5, 2.7, 0.8);
    Joint3D rf(pt_rf, un);
    Skeleton3D rightfront(rf, vide, "rf");
    rightfront.add_child(xrightfront);
    Point3f pt_ne(3.8, 0, 0);
    Joint3D ne(pt_ne, un);
    Skeleton3D neck(ne, vide, "neck");
    neck.add_constraint(0, 1);                     // keep the back horizontal
    neck.add_child(leftfront); neck.add_child(rightfront); neck.add_child(head);
    
    // back legs
    Point3f pt_xxlb(0.1, 1, 0);
    Joint3D xxlb(pt_xxlb, un);
    Skeleton3D xxleftback(xxlb, vide, "xxlb");
    Point3f pt_xxrb(-0.1, 1, 0);
    Joint3D xxrb(pt_xxrb, un);
    Skeleton3D xxrightback(xxrb, vide, "xxrb");
    Point3f pt_xlb(0.5, 1.4, 0);
    Joint3D xlb(pt_xlb, un);
    Skeleton3D xleftback(xlb, vide, "xlb");
    xleftback.add_child(xxleftback);
    Point3f pt_xrb(-0.5, 1.4, 0);
    Joint3D xrb(pt_xrb, un);
    Skeleton3D xrightback(xrb, vide, "xrb");
    xrightback.add_child(xxrightback);
    Point3f pt_lb(.5, 2.7, -0.8);
    Joint3D lb(pt_lb, un);
    Skeleton3D leftback(lb, vide, "lb");
    leftback.add_child(xleftback);
    Point3f pt_rb(-.5, 2.7, 0.8);
    Joint3D rb(pt_rb, un);
    Skeleton3D rightback(rb, vide, "rb");
    rightback.add_child(xrightback);

    // back legs
    /*Point3f pt_xxblf(0.1, 1, 0);
    Joint3D xxblf(pt_xxblf, un);
    Skeleton3D xxbleft(xxblf, vide, "xxblf");
    Point3f pt_xxbrf(-0.1, 1, 0);
    Joint3D xxbrf(pt_xxbrf, un);
    Skeleton3D xxbright(xxbrf, vide, "xxbrf");
    Point3f pt_xlbf(-0.3, 1, 0);
    Joint3D xlbf(pt_xlbf, un);
    Skeleton3D xleftback(xlbf, vide, "xlb");
    xleftback.add_child(xxbleft);
    Point3f pt_xrbf(-1, 1, 0);
    Joint3D xrbf(pt_xrbf, un);
    Skeleton3D xrightback(xrbf, vide, "xrb");
    xrightback.add_child(xxbright);
    Point3f pt_lbf(0, 2.88, -1);
    Joint3D lbf(pt_lbf, un);
    Skeleton3D leftback(lbf, vide, "lb");
    leftback.add_child(xleftback);
    Point3f pt_rbf(-1, 2.7, 1);
    Joint3D rbf(pt_rbf, un);
    Skeleton3D rightback(rbf, vide, "rb");
    rightback.add_child(xrightback);*/
    
    // tail
    Point3f pt_xt(-0.3, 1.2, 0);
    Joint3D xt(pt_xt, un);
    Skeleton3D xtail(xt, vide, "xtail");
    Point3f pt_t(-0.5, 1.8, 0);
    Joint3D t(pt_t, un);
    Skeleton3D tail(t, vide, "tail");
    tail.add_child(xtail);
    
    Point3f pt_hip_p(0, 0, 0);
    Joint3D hip_p(pt_hip_p, un);
    Skeleton3D hip(hip_p, vide, "hip");
    hip.add_child(neck);hip.add_child(leftback);hip.add_child(rightback);hip.add_child(tail);
    return hip;
}


Skeleton3D test_bird(){
    vector<Skeleton3D> vide;
    Quaternion un(1, 0, 0, 0);
    Point3f pt_xhead(-0.3, -0.3, 0);
    Joint3D j_xhead(pt_xhead, un);
    Skeleton3D xhead(j_xhead, vide, "xhead");
    Point3f pt_head(-0.5, -0.5, 0);
    Joint3D j_head(pt_head, un);
    Skeleton3D head(j_head, vide, "head");
    head.add_child(xhead);
    Point3f pt_xt(0.4, 0.4, 0);
    Joint3D j_xt(pt_xt, un);
    Skeleton3D xtail(j_xt, vide, "xtail");
    Point3f pt_t(0.5, 0.5, 0);
    Joint3D j_t(pt_t, un);
    Skeleton3D tail(j_t, vide, "tail");
    tail.add_child(xtail);
    
    Point3f pt_xl(-0.7, 0.5, 1);
    Joint3D j_xl(pt_xl, un);
    Skeleton3D xl(j_xl, vide, "extreme_left");
    Point3f pt_l(-0.7, 0.5, 1);
    Joint3D j_l(pt_l, un);
    Skeleton3D l(j_l, vide, "left");
    l.add_child(xl);
    
    Point3f pt_xr(-0.5, 0.7, -1);
    Joint3D j_xr(pt_xr, un);
    Skeleton3D xr(j_xr, vide, "extreme_right");
    Point3f pt_r(-0.5, 0.7, -1);
    Joint3D j_r(pt_r, un);
    Skeleton3D r(j_r, vide, "right");
    r.add_child(xr);
    
    Point3f pt_hip(1, 1, 1);
    Joint3D j_hip(pt_hip, un);
    Skeleton3D hip(j_hip, vide, "hip");
    hip.add_child(l);hip.add_child(r);hip.add_child(head);hip.add_child(tail);
    return hip;
}
