//
//  Skeleton2D.cpp
//  projetc
//
//  Created by Coralie DAVID on 27/10/2017.
//  Copyright © 2017 Manon Romain. All rights reserved.
//

#include "Skeleton2D.hpp"
#define PI 3.14159265

void Skeleton2D::transform(Mat& h){
    //cout<<"transformation de "<<name<<" : "<<this->pos;
    Point3f p = h*(this->pos);
    this->setPos(p);
    //cout<<" devient "<<this->pos<<endl;
    for (size_t i=0; i<children.size(); i++){
        children[i].transform(h);
    }
    updateMinMax();
}
void Skeleton2D::updateMinMax(){
    pos = pos/pos.z;
    double minx, maxx, miny, maxy;
    minx = pos.x; maxx = pos.x;
    miny = pos.y; maxy = pos.y;
    for (size_t i=0; i<children.size(); i++){
        minx = min(minx, children[i].minMax[0]);
        maxx = max(maxx, children[i].minMax[1]);
        miny = min(miny, children[i].minMax[2]);
        maxy = max(maxy, children[i].minMax[3]);
    }
    this->minMax[0] = minx;
    this->minMax[1] = maxx;
    this->minMax[2] = miny;
    this->minMax[3] = maxy;
}
void Skeleton2D::normalize(int width, int height){
    Mat m = Mat::zeros(3, 3, CV_32FC1);
    //cout<<s<<endl;
    
    //cout<<"Min et max :"<<minMax[0]<<", "<<minMax[1]<<", "<<minMax[2]<<", "<<minMax[3]<<endl;
    float scale = min(width, height)/(max(minMax[1]-minMax[0],minMax[3]-minMax[2])+4);
    m.at<float>(0,0) = scale;
    m.at<float>(1,1) = scale;
    m.at<float>(2,2) = 1;
    m.at<float>(0,2) = scale*(-minMax[0]);
    m.at<float>(1,2) = scale*(-minMax[2]);
    //cout<<"Matrice de normalisation : \n"<<m<<endl;
    transform(m);
}

Mat Skeleton2D::toMat(int width, int height){
    Mat m = Mat::zeros(width, height, CV_8U);
    Mat debug_im = 255*Mat::ones(width, height, CV_8U);
    auxMat(width, height, m, debug_im);
    //imshow("skeleton", debug_im); waitKey(0);
    return m;
};

void Skeleton2D::auxMat(int width, int height, Mat& dest, Mat& debug_im){
    Point3f p = this->pos;
    Size pt(2, 2);
    Point p2(p.x/p.z, p.y/p.z);
    //ellipse(dest, p2, pt, 0, 0, 360, Scalar(255, 0, 255), -1);
    for (size_t i=0; i<children.size(); i++){
        oneMat(width, height, dest, debug_im, i);
        children[i].auxMat(width, height, dest, debug_im);
    }
};

void Skeleton2D::oneMat(int width, int height, Mat& dest, Mat& debug_im, int i){
    Point3f p = this->pos;
    Size pt(2, 2);
    Point p2(p.x/p.z, p.y/p.z);
    Point p_dest(children[i].pos.x/children[i].pos.z, children[i].pos.y/children[i].pos.z);
    Point barycenter = (p2+p_dest)/2; Point diff = (p2-p_dest);
    
    Point center(barycenter.x, barycenter.y);
    Point vector(diff.x, diff.y);
    //cout<<name<<" to "<<children[i].name<<endl;
    double d = norm(diff);
    double theta;
    
    if (vector.x==0){
        theta = PI/2;
        //cout<<" x nul"<<endl;
    }
    
    else {
        theta = atan((float)vector.y/vector.x);
    }
    if (children[i].name=="neck"){
        Size s(d/2, 3./5*d/2);
        //Size s(10, 10);
        Point pt1 = p2;
        Point pt2 = p_dest + Point(0, 5. / 6 * d); // + Point(100,100);
        //cout << "point: " << d << pt1 << pt2 << endl;
        rectangle(dest, pt1, pt2, Scalar(255, 255, 255), cv::FILLED);
        //ellipse(dest, center, s, theta*180/PI, 0, 360, Scalar(255, 0, 0), -1);
    }
    else {
        //Size s(d/2, pow(d, 3./2)/50);
        Size s(d / 2, 10);
        ellipse(dest, center, s, theta*180/PI, 0, 360, Scalar(255, 0, 0), cv::FILLED);
    }
    //ellipse(debug_im, center, pt, 0, 0, 360, Scalar(0, 0, 0), -1);
    arrowedLine(debug_im, p2, p_dest, Scalar(0, 0, 0));
    putText(debug_im, children[i].name, center, FONT_HERSHEY_SIMPLEX, 0.4, Scalar(0,0, 0));
};


