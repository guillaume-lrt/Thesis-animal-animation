//
//  shape_context.hpp
//  projetc
//
//  Created by Coralie DAVID on 05/12/2017.
//  Copyright © 2017 Manon Romain. All rights reserved.
//

#ifndef shape_context_hpp
#define shape_context_hpp

#include <stdio.h>
#include <vector>
#include <stack>
#include <iostream>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/shape/shape.hpp>
#include <opencv2/videoio/videoio.hpp>
#include <opencv2/highgui/highgui.hpp>


using namespace std;
using namespace cv;

class ShapeContextDistance{
    Ptr<ShapeContextDistanceExtractor> al;
    public :
    ShapeContextDistance(){
        al = createShapeContextDistanceExtractor();
        al->setShapeContextWeight(1);
        al->setBendingEnergyWeight(0);
        al->setImageAppearanceWeight(0);
    }
    double distance(vector<Point> c1, vector<Point> c2);
    vector<vector<double>> apply(vector<Mat> vm);

};

#endif /* shape_context_hpp */
