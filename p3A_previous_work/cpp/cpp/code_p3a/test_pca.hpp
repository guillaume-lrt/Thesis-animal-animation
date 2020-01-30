//
//  test_pca.hpp
//  projetc
//
//  Created by Coralie DAVID on 21/11/2017.
//  Copyright © 2017 Manon Romain. All rights reserved.
//

#ifndef test_pca_hpp
#define test_pca_hpp

#include <stdio.h>
#include <iostream>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/videoio/videoio.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;
using namespace cv;

vector<Point> testContours(Mat& m);
double get_angle_pca(Mat& m);
Mat get_rot(Mat& m1, Mat& m2);
Mat get_trans(Mat& m1, Mat& m2);
void drawAxis(Mat& img, Point p, Point q, Scalar colour, const float scale = 0.2);

#endif /* test_pca_hpp */
