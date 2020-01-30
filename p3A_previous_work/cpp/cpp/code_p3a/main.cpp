//
//  main.cpp
//  projetc
//
//  Created by Coralie DAVID on 27/10/2017.
//  Copyright © 2017 Manon Romain. All rights reserved.
//

#include <iostream>
#include <filesystem>
#include <fstream>
#include <stdio.h>
//#include <unistd.h>
#include <queue>

#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/videoio/videoio.hpp>
#include <opencv2/highgui/highgui.hpp>


#include "Data_Structure/Quaternion.hpp"
#include "Data_Structure/Joint3D.hpp"
#include "Data_Structure/Skeleton3D.hpp"
#include "Data_Structure/Skeleton2D.hpp"
#include "utils.hpp"
#include "test_pca.hpp"
#include "test_skeletons.hpp"

#include "animating.hpp"
#include "Data_Structure/Sprite.hpp"
using namespace cv;
//using namespace std;
#include "string"


Mat optimize(Mat& im, Skeleton2D& t){
    Mat c, c2, c3;
    c = t.toMat(im.rows, im.cols);
    
    //Min Max Scaling using bboxes
    Mat h_rect = find_rect_match(im, c);
    t.transform(h_rect);
    c2 = t.toMat(im.rows, im.cols);
    
    //Getting rot and trans using pca
    //Mat h = get_rot(im, c2);
    Mat h = get_trans(im, c2);
    t.transform(h);
    c3 = t.toMat(im.rows, im.cols);
    
    return c3;
}
void show_merged(const String name, const Mat& m1, const Mat& m2, String strr="None"){
    Mat i = 255*Mat::ones(m1.rows, m1.cols, m1.type());
    vector<Mat> mix = {i-m1/4-m2/2, i-m1/4, i-m1/2-m2/2};
    Mat merging, out;
    merge(mix, merging);
    
    
    resize(merging, out, Size(300, 300));
    //out = merging;
    Point center(0, 25);
    putText(out, "r="+strr, center, FONT_HERSHEY_SIMPLEX, 1, Scalar(0, 0, 255));
    imshow(name, out);
    return;
}
double try_and_optimize(Mat& target, Vec3f& v, Skeleton3D& s, double angle_min, double angle_max){
    Mat c;

    int n = 10;
    double var_k = (angle_max-angle_min)/ (double) n;
    
    queue<Skeleton3D*> to_do;
    to_do.push(&s);
    
    v = v/norm(v);
    Mat rep = repere_plane(v);
    double r_max = 0;
    //Breadth First Search
    while (!to_do.empty()){
        Skeleton3D* moving = to_do.front();
        to_do.pop();
        for (int j=0; j<moving->get_children_size(); j++){
            to_do.push(moving->get_child(j));
        }
        
        double k_max = 0;
        double angle = angle_min;
        r_max = 0;
        for (int k=0; k<=n; k++){
            Quaternion q(angle, v);
            moving->rotate(q);
            s.updateAbsolutePosition();
            
            Skeleton2D t = s.project(rep);
            t.normalize(target.rows, target.cols);
            c = optimize(target, t);
            //imshow("c", c); waitKey();
            Scalar r = iou(target, c);
            if (r[0]>r_max){
                k_max = k;
                r_max = r[0];
            }
            //continue
            angle = var_k;
        }
        //back to optimal angle
        Quaternion q((k_max-n)*var_k, v);
        moving->rotate(q);
        s.updateAbsolutePosition();
        
    }
    cout<<r_max<<endl;
    show_merged("m2", target, c);
    return r_max;
}

double new_conf(Mat& target, Vec3f& v, Mat& rep, Skeleton3D& s, double angle_max){
    Mat c;
    Scalar r;
    
    //symetry
    /*double alea = (double) rand()/RAND_MAX;
    if (alea<0.01){
        cout<<"yes"<<endl;
        Mat h = Mat::zeros(3, 3, CV_32F);
        h.at<float>(0, 0) = 1; h.at<float>(1, 1) = -1; h.at<float>(2, 2) = 1;
        s.transform(h);
    }*/
    
    
    queue<Skeleton3D*> to_do;
    to_do.push(&s);
    //bool first = true;
    //Breadth First
    while (!to_do.empty()){
        Skeleton3D* moving = to_do.front();
        to_do.pop();
        for (int j=0; j<moving->get_children_size(); j++){
            to_do.push(moving->get_child(j));
        }
        double angle = 2*angle_max*PI/180*(double)rand() / RAND_MAX - angle_max*PI/180;
        //if (!first){
        //cout << "Random angle: " << angle << moving->get_name() << endl;
        //double angle = 0.137607;
        Quaternion q(angle, v);
        moving->rotate(q);
        //}
        //else
        //    first = false;
        s.updateAbsolutePosition();
        
        Skeleton2D t = s.project(rep);
        t.normalize(target.rows, target.cols);
        c = optimize(target, t);
        r = iou(target, c);
    }
    
    return r[0];
}

void recuit_simule(Mat& target, Vec3f& v, Skeleton3D& s_max){
    // target = mask/1.png

    double T = 20;
    double decr = 0.999;
    double limit = 0.0001;
    double angle_max = 8;       // = 8;
    double k_b = 2*1e-4;
    
    v = v/norm(v);
    Mat rep = repere_plane(v);
    
    char key = 'a';
    ofstream fichier("energy.txt", ios::out | ios::trunc);  //déclaration du flux et ouverture du fichier
    
    if(fichier)  // si l'ouverture a réussi
    {
        double r, r_max=0;
        int j = 0;
        while (T>limit && key!='q'){
            int new_confs = 0;
            cout << "New conf: " << new_confs << ", and j: " << j << endl;
            //int new_confs_alea = 0;
            for (int i=0; i<10; i++){
                Skeleton3D s(s_max);
                r = new_conf(target, v, rep, s, angle_max);
                if (j==0){
                    k_b = r/(100*T);
                }
                //double alea = (double) rand()/RAND_MAX;
                if (r>r_max) {
                    fichier<<"0;";
                    r_max = r;
                    s_max = s;
                }
//                else if (alea<exp(-(r_max-r)/(k_b*T))){
//                    fichier<<"1;";
//                    new_confs_alea++;
//                    r_max = r;
//                    s_max = s;
//                }
                else {
                    fichier<<"2;";
                    new_confs++;
                }
                fichier<<r_max<<";"<<T<<endl;
            }
            //cout<<"ok:"<<new_confs_alea*100/(new_confs_alea+new_confs)<<"%"<<endl;
            if (j%2==0){
                //cout<<"T="<<T<<endl;
                
                Skeleton2D t = s_max.project(rep);
                t.normalize(target.rows, target.cols);
                Mat c = optimize(target, t);
                show_merged("recuit", target, c, to_string(r_max));//+to_string(j)
                key = waitKey(50);
            }
            T *= decr;
            j++;
        }
        fichier.close();  // on referme le fichier
    }
    else  // sinon
        cerr << "Erreur à l'ouverture !" << endl;
}


void par_morceaux(Mat& target, Vec3f& v, Skeleton3D& s_max){

    v = v/norm(v);
    Mat rep = repere_plane(v);
    
    char key = 'a';
    
    
    
}

int main(int argc, const char * argv[]) {
    
    //get cwd
    /*size_t s;
    char* buffer;
    buffer = getcwd(buffer, s);
    cout<<buffer<<endl;*/
    
    //Open and display target im
    int n = 7;
    std::string path_temp = argv[0];
    String file_path = path_temp + "/../../../cpp/data/elephants2_mask/";
    //String file_path = path_temp + "/../../../cpp/data/original/";
    //String file_path = "../../../data/original/";
    //String file_path = "../../../data/birds_masks/";
    //String file_path = "../../../data/oies_masks/";
    //String file_path = "../../../data/oies_bas_masks/";
    //String file_path = "../../../data/oies2_masks/";
    SnapshotGraph g(file_path, n, true);

    vector<int> path = g.get_path();
    Vec3f v_i(1, 0, 5);
    v_i = v_i/norm(v_i);
    Sprite s(file_path, path, false, v_i);
    char key = 'a';
    int theta = 30;
    while (key!='q'){
        s.display_one_horiz((double)theta);
        s.display_one_vert((double)theta); key = waitKey();
        if ((int)key==0)
            theta = (theta+5)%360;
        if ((int)key==1)
            theta = (theta-5)%360;
    }
    
    int i_query = 1;
    stringstream query;
    query << file_path << i_query << ".png";
    Mat im = imread(query.str(), cv::IMREAD_GRAYSCALE);
    
    resize(im, im, Size(), 1.6, 1.6);
    //imshow("hello", im);waitKey(0);
    cout<<"Image: " << query.str() << ", with shape = "<<im.size<<endl;
    
    //Test 3D skeleton
    Skeleton3D hip = test3D();
    hip.updateAbsolutePosition();
    Vec3f v(0, 0, 1);
    v = v/norm(v);
    Mat rep = repere_plane(v);

    cout << "debug1" << endl;

    recuit_simule(im, v, hip);

    cout << "Key: " << key << endl;

    //double r_prec = try_and_optimize(im, v, hip, -.3, .3); waitKey();
    //double r = try_and_optimize(im, v, hip, -.3, .3); waitKey();
    //while (r_prec!=r){
    //    r_prec = r;
    //    r = try_and_optimize(im, v, hip, -.3, .3); waitKey();
    //}
    im.release();
    return 0;
}
