//
//  shape_context.cpp
//  projetc
//
//  Created by Coralie DAVID on 05/12/2017.
//  Copyright © 2017 Manon Romain. All rights reserved.
//

#include "shape_context.hpp"
vector<Point> simpleContour(const Mat& currentQuery, int n=400)
{
    vector<vector<Point> > _contoursQuery;
    vector <Point> contoursQuery;
    findContours(currentQuery, _contoursQuery, RETR_LIST, CHAIN_APPROX_NONE);
    for (size_t border=0; border<_contoursQuery.size(); border++)
    {
        for (size_t p=0; p<_contoursQuery[border].size(); p++)
        {
            contoursQuery.push_back( _contoursQuery[border][p] );
        }
    }
    
    // In case actual number of points is less than n
    int dummy=0;
    if ((int)contoursQuery.size()<n)
        cout<<"Adding "<<n-contoursQuery.size()<<" points"<<endl;
    for (int add=(int)contoursQuery.size()-1; add<n; add++)
    {
        //cout<<add<<endl;
        contoursQuery.push_back(contoursQuery[dummy++]); //adding dummy values
    }
    
    // Uniformly sampling
    random_shuffle(contoursQuery.begin(), contoursQuery.end());
    vector<Point> cont;
    
    bool debug = true;
    Mat m;
    if (debug){
        cvtColor(currentQuery, m, cv::COLOR_GRAY2BGR);
    }
    int i = 0;
    double criterium = 2;
    while (cont.size()<n)
    {
        if (i>=contoursQuery.size()) {
            criterium *= 0.8;
        }
        i = i % contoursQuery.size();
        if (debug && cont.size()%200==0) {imshow("opencv contour sampling", m); waitKey(0);}
        bool ok = true;
        for (Point p: cont)
            if (norm(p-contoursQuery[i])<criterium){
                ok = false; break;
            }
        if (ok){
            cont.push_back(contoursQuery[i]);
            if (debug){
                circle(m, contoursQuery[i], 1, Scalar(0, 0, 255), -1);
            }
            
        }
        else {
            if (debug){
                circle(m, contoursQuery[i], 1, Scalar(0, 255, 0), -1);
            }
        }
        
        i++;
    }
    if (debug)
        imshow("opencv contour sampling", m); waitKey();

    return cont;
}


double ShapeContextDistance::distance(vector<Point> contour1, vector<Point> contour2){
    double d = this->al->computeDistance(contour1, contour2);
    return d;
}

vector<vector<double>> ShapeContextDistance::apply(vector<Mat> vm){
    int n = vm.size();
    vector<vector<double>> g(n);
    for (int i=0; i<n; i++){
        vector<double> gi(n);
        g[i] = gi;
        g[i][i] = 0;
    }
    vector<vector<Point>> contours;
    int i=0;
    for (Mat m:vm){
        contours.push_back(simpleContour(m));
        for (int j=0; j<i; ++j){
            double s = distance(contours[i], contours[j]);
            cout<<"Distance entre "<<i<<" et "<<j<<" = "<<s<<endl;
            //cout<<s<<",";
            g[i][j] = s; g[j][i] = s;
        }
        i++;
    }
    return g;
}
