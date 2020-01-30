//
//  animating.hpp
//  projetc
//
//  Created by Coralie DAVID on 17/12/2017.
//  Copyright © 2017 Manon Romain. All rights reserved.
//

#ifndef animating_hpp
#define animating_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/shape/shape.hpp>
#include <opencv2/videoio/videoio.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "shape_context.hpp"
#include "constantes.hpp"
#include "utils.hpp"


using namespace std;
using namespace cv;

class SnapshotGraph{
    
    double lambda_1 = 2.5;
    double lambda_2 = 0.5;
    bool full_cycle;
    
    bool verbose = true;
    
    vector<vector<double>> g;
    double alpha = 1;
    int start, end;
    vector<int> path;
    vector<bool> outliers;
    private:
        double distance(int i, int j);
        vector<int> generate_path(int l);
        vector<int> generate_path_half_cycle(int l);
        vector<int> generate_path_full_cycle(int l);
        double energy_path(int l, vector<int> path);
        void load_from_txt(String file);
        void write_in_file(String file);
    public:
        SnapshotGraph(String path, int n, bool cycle);
        void delete_outliers(double c);
        void select_min_max();
        void create_path();
        vector<int> get_path(){
            return this->path;
        }
    
};




#endif /* animating_hpp */
