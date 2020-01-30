//
//  utils.h
//  projetc
//
//  Created by Coralie DAVID on 04/11/2017.
//  Copyright © 2017 Manon Romain. All rights reserved.
//

#ifndef utils_h
#define utils_h
#include <stdio.h>
#include <iostream>
#include <vector>
using namespace std;

#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/videoio/videoio.hpp>
#include <opencv2/highgui/highgui.hpp>
using namespace cv;

Point3f operator*(Mat& M, const Point3f& p);

Scalar iou(const Mat& M1, const Mat& M2);
Mat find_rect_match(Mat& m1, Mat& m2);
Mat transform_rect(Mat& m1, Mat& m2);
Mat transform_im(Mat& h, Mat& m2);
Mat translation_matrix(double dx, double dy);
Mat rotateY(double theta);
Mat repere_plane(const Vec3f& d);
Mat repere_plane23(const Vec3f& d);

double avg (vector<double>& v);
double avg (vector<vector<double>>& v);

double variance (vector<double>& v , double mean);
double variance (vector<vector<double>>& v , double mean);

const vector<String> explode(const String& s, const char& c);


#endif /* utils_h */
