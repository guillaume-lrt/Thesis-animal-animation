//
//  test_pca.cpp
//  projetc
//
//  Created by Coralie DAVID on 21/11/2017.
//  Copyright © 2017 Manon Romain. All rights reserved.
//

#include "test_pca.hpp"
#include <opencv2/core/core_c.h>

RNG rng(12345);
bool DEBUGGING = false;

vector<Point> testContours(Mat& src){
    int thresh = 10;
    Mat src_gray, canny_output;
    vector<vector<Point> > contours;
    vector<Point> contours_flat;
    vector<Vec4i> hierarchy;
    if (src.channels()==3){
        cvtColor(src, src_gray, cv::COLOR_BGR2GRAY);
    }
    else {
        src_gray = src;
    }
        
    /// Detect edges using canny
    Canny( src_gray, canny_output, thresh, thresh*2, 3 );
    /// Find contours
    findContours( canny_output, contours, hierarchy, cv::RETR_LIST, cv::CHAIN_APPROX_SIMPLE, Point(0, 0) );
    
    /// Draw contours
    Mat drawing = Mat::zeros( canny_output.size(), CV_8UC3 );
    for( int i = 0; i< contours.size(); i++ )
    {
        //cout<<"Contour n°"<<i<<" de taille "<<contours[i].size()<<endl;
        if (contours[i].size()>10){
            //Scalar color = Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) );
            //drawContours( drawing, contours, i, color, 2, 8, hierarchy, 0, Point() );
            for( int j = 0; j< contours[i].size(); j++ )
                contours_flat.push_back(contours[i][j]);
        }
        
    }
    
    /// Show in a window
    //namedWindow( "Contours", CV_WINDOW_AUTOSIZE );
    //imshow( "Contours", drawing );
    return contours_flat;
}

double get_angle_pca(Mat& m){
    vector<Point> pts = testContours(m);
    //Construct a buffer used by the pca analysis
    int sz = static_cast<int>(pts.size());
    Mat data_pts = Mat(sz, 2, CV_64FC1);
    for (int i = 0; i < data_pts.rows; ++i)
    {
        data_pts.at<double>(i, 0) = pts[i].x;
        data_pts.at<double>(i, 1) = pts[i].y;
    }
    //Perform PCA analysis
    PCA pca_analysis(data_pts, Mat(), CV_PCA_DATA_AS_ROW);
    //Store the center of the object
    Point cntr = Point(static_cast<int>(pca_analysis.mean.at<double>(0, 0)),
                       static_cast<int>(pca_analysis.mean.at<double>(0, 1)));
    //Store the eigenvalues and eigenvectors
    
    vector<Point2d> eigen_vecs(2);
    vector<double> eigen_val(2);
    for (int i = 0; i < 2; ++i)
    {
        eigen_vecs[i] = Point2d(pca_analysis.eigenvectors.at<double>(i, 0),
                                pca_analysis.eigenvectors.at<double>(i, 1));
        eigen_val[i] = pca_analysis.eigenvalues.at<double>(0, i);
    }
    if (DEBUGGING){
        // Draw the principal components
        circle(m, cntr, 3, Scalar(255, 0, 255), 2);
        Point p1 = cntr + 0.002 * Point(static_cast<int>(eigen_vecs[0].x * eigen_val[0]), static_cast<int>(eigen_vecs[0].y * eigen_val[0]));
        Point p2 = cntr - 0.002 * Point(static_cast<int>(eigen_vecs[1].x * eigen_val[1]), static_cast<int>(eigen_vecs[1].y * eigen_val[1]));
        drawAxis(m, cntr, p1, Scalar(0, 255, 0), 1);
        drawAxis(m, cntr, p2, Scalar(255, 255, 0), 5);
    }
    return atan2(eigen_vecs[0].y, eigen_vecs[0].x);;
    
}

Mat get_rot(Mat& m1, Mat& m2){

    double theta1 = get_angle_pca(m1);
    double theta2 = get_angle_pca(m2);
    double theta = theta1 - theta2;
    Mat h = Mat::zeros(3, 3, CV_32FC1);
    h.at<float>(0,0) = cos(theta);
    h.at<float>(1,1) = cos(theta);
    h.at<float>(2,2) = 1;
    h.at<float>(0,1) = -sin(theta);
    h.at<float>(1,0) = sin(theta);
    Moments mom1 = moments(m1);
    Moments mom2 = moments(m2);
    h.at<float>(0,2) = (mom1.m10/mom1.m00 - cos(theta)*mom2.m10/mom2.m00 + sin(theta)*mom2.m01/mom2.m00);
    h.at<float>(1,2) = (mom1.m01/mom1.m00 - sin(theta)*mom2.m10/mom2.m00 - cos(theta)*mom2.m01/mom2.m00);
    return h;
}

Mat get_trans(Mat& m1, Mat& m2){
    Mat h = Mat::zeros(3, 3, CV_32FC1);
    h.at<float>(0,0) = 1;
    h.at<float>(1,1) = 1;
    h.at<float>(2,2) = 1;
    Moments mom1 = moments(m1);
    Moments mom2 = moments(m2);
    h.at<float>(0,2) = (mom1.m10/mom1.m00 - mom2.m10/mom2.m00);
    h.at<float>(1,2) = (mom1.m01/mom1.m00 - mom2.m01/mom2.m00);
    return h;
}


void drawAxis(Mat& img, Point p, Point q, Scalar colour, const float scale)
{
    double angle;
    double hypotenuse;
    angle = atan2( (double) p.y - q.y, (double) p.x - q.x ); // angle in radians
    hypotenuse = sqrt( (double) (p.y - q.y) * (p.y - q.y) + (p.x - q.x) * (p.x - q.x));
    //    double degrees = angle * 180 / CV_PI; // convert radians to degrees (0-180 range)
    //    cout << "Degrees: " << abs(degrees - 180) << endl; // angle in 0-360 degrees range
    // Here we lengthen the arrow by a factor of scale
    q.x = (int) (p.x - scale * hypotenuse * cos(angle));
    q.y = (int) (p.y - scale * hypotenuse * sin(angle));
    line(img, p, q, colour, 1, cv::LINE_AA);
    // create the arrow hooks
    p.x = (int) (q.x + 9 * cos(angle + CV_PI / 4));
    p.y = (int) (q.y + 9 * sin(angle + CV_PI / 4));
    line(img, p, q, colour, 1, cv::LINE_AA);
    p.x = (int) (q.x + 9 * cos(angle - CV_PI / 4));
    p.y = (int) (q.y + 9 * sin(angle - CV_PI / 4));
    line(img, p, q, colour, 1, cv::LINE_AA);
}
