//
//  utils.cpp
//  projetc
//
//  Created by Coralie DAVID on 03/11/2017.
//  Copyright © 2017 Manon Romain. All rights reserved.
//


#include "utils.hpp"

Point3f operator*(Mat& M, const Point3f& p)
{
    Mat_<float> src(3/*rows*/,1 /* cols */);
    
    src(0,0)=p.x;
    src(1,0)=p.y;
    src(2,0)=p.z;
    
    Mat_<float> dst = M*src; //USE MATRIX ALGEBRA
    return Point3f(dst(0,0),dst(1,0), dst(2,0));
}

Scalar iou(const Mat& M1, const Mat& M2){
    M1.convertTo(M1, CV_8U);
    M2.convertTo(M2, CV_8U);
    Mat i, u;
    bitwise_and(M1, M2, i);
    bitwise_or(M1, M2, u);
    Scalar r = sum(i)/sum(u);
    return r;
};

Mat find_rect_match(Mat& m1, Mat& m2){
    Mat h = Mat::zeros(3, 3, CV_32FC1);
    Rect r1 = boundingRect(m1);
    Rect r2 = boundingRect(m2);
    float s = min(((float)r1.width)/r2.width, ((float)r1.height)/r2.height);
    h.at<float>(0,0) = s;
    h.at<float>(1,1) = s;
    h.at<float>(0,2) = s*(r1.x-r2.x);
    h.at<float>(1,2) = s*(r1.y-r2.y);
    h.at<float>(2,2) = 1;
    return h;
};

Mat transform_rect(Mat& m1, Mat& m2){
    Mat h = find_rect_match(m1, m2);
    Mat c2;
    warpPerspective(m2, c2, h, c2.size());
    return c2;
}

Mat transform_im(Mat& h, Mat& m2){
    cout<<"MATRICE de rotation :"<<h<<endl;
    Mat c2;
    warpPerspective(m2, c2, h, c2.size());
    return c2;
}

Mat repere_plane(const Vec3f& d){
    Mat m = Mat::zeros(3, 3, CV_32FC1);
    Vec3f i,j;
    if (d[0]==0)
        i = Vec3f(1, 0, 0);
    else
        i = Vec3f(-d[1]/d[0], 1, 0);
    i = i/norm(i);
    j = d.cross(i);
    m.at<float>(0, 0) = i[0]; m.at<float>(1, 0) = j[0]/norm(j); m.at<float>(2, 0) = d[0];
    m.at<float>(0, 1) = i[1]; m.at<float>(1, 1) = j[1]/norm(j); m.at<float>(2, 1) = d[1];
    m.at<float>(0, 2) = i[2]; m.at<float>(1, 2) = j[2]/norm(j); m.at<float>(2, 2) = d[2];
    return m;
}

Mat repere_plane23(const Vec3f& d){
    Mat m = Mat::zeros(2, 3, CV_32FC1);
    Vec3f i,j;
    if (d[0]==0)
        i = Vec3f(1, 0, 0);
    else
        i = Vec3f(-d[1]/d[0], 1, 0);
    i = i/norm(i);
    j = d.cross(i);
    m.at<float>(0, 0) = i[0]; m.at<float>(1, 0) = j[0]/norm(j);
    m.at<float>(0, 1) = i[1]; m.at<float>(1, 1) = j[1]/norm(j);
    m.at<float>(0, 2) = i[2]; m.at<float>(1, 2) = j[2]/norm(j);
    return m;
}

Mat rotateY(double theta){
    Mat m = Mat::zeros(3, 3, CV_32FC1);
    m.at<float>(0, 0) = cos(theta); /*******0;***********/  m.at<float>(2, 0) = -sin(theta);
    /*************0;**************/ m.at<float>(1, 1) = 1;/***************0;**************/
    m.at<float>(2, 0) = -sin(theta); /*******0;***********/ m.at<float>(2, 2) = cos(theta);
    return m;
}

Mat translation_matrix(double dx, double dy){
    Mat m = Mat::zeros(3, 3, CV_32FC1);
    m.at<float>(0, 0) = 1; m.at<float>(0, 2) = dx;
    m.at<float>(1, 1) = 1; m.at<float>(1, 2) = dy;
    m.at<float>(2, 2) = 1;
    return m;
}

//Function for average
double avg ( vector<double>& v )
{
    double return_value = 0.0;
    double n = v.size();
    
    for ( int i=0; i < n; i++)
    {
        return_value += v[i];
    }
    //cout<<"r="<<return_value/n;
    
    return return_value/n;
}

double avg ( vector<vector<double>>& v )
{
    double return_value = 0.0;
    int n = 0;
    
    for ( int i=0; i < v.size(); i++)
    {
        int ni = v[i].size();
        double a = avg(v[i]);
        //cout<<" = "<<a*ni<<endl;
        return_value += ni*a;
        n += ni;
    }
    
    return ( return_value / n);
}
//****************End of average funtion****************


//Function for variance
double variance ( vector<double>& v , double mean )
{
    double sum = 0.0;
    double temp = 0.0;
    
    for ( int j =0; j < v.size(); j++)
    {
        temp = pow((v[j] - mean),2);
        sum += temp;
    }
    
    return sum/(double)v.size();
}

double variance (vector<vector<double>>& v , double mean)
{
    double sum = 0.0;
    int n = 0;
    
    for ( int j =0; j < v.size(); j++)
    {
        n += v[j].size();
        for (double a: v[j])
            sum += pow((a - mean),2);
    }
    
    return sum/n;
}

const vector<String> explode(const String& s, const char& c)
{
    String buff{""};
    vector<String> v;
    
    for(auto n:s)
    {
        if(n != c) buff+=n; else
            if(n == c && buff != "") { v.push_back(buff); buff = ""; }
    }
    if(buff != "") v.push_back(buff);
    
    return v;
}
