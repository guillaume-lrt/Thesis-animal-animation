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
#include <math.h>

#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/videoio/videoio.hpp>
#include <opencv2/highgui/highgui.hpp>


#include "Data_Structure/Quaternion.hpp"
#include "Data_Structure/Joint3D.hpp"
#include "Data_Structure/Skeleton3D.hpp"
#include "Data_Structure/Skeleton2D.hpp"
#include "Data_Structure/utils.hpp"
#include "test_pca.hpp"
#include "test_skeletons.hpp"

#include "animating.hpp"
#include "Data_Structure/Sprite.hpp"
using namespace cv;
//using namespace std;
#include "string"
#include <map>

Mat rotation_matrix(Vec3f n, double theta) {
    // n: rotation axis; =(n_x,n_y,n_z)
    // theta: angle, perpendicular to the axis
    Mat res(Size(3, 3), CV_32FC1);

    res.at<float>(0, 0) = cos(theta) + n[0]*n[0]*(1-cos(theta));
    res.at<float>(0, 1) = n[0]*n[1]*(1-cos(theta)) - n[3]*sin(theta);
    //...

    return res;
}

Mat rotation_matrix(double theta) {
    // We assume that all rotations are along the z axis, i.e n = (0,0,1)
    Mat res(Size(3, 3), CV_32FC1);
    res.at<double>(0, 0) = cos(theta);
    res.at<double>(0, 1) = -sin(theta);
    res.at<double>(0, 2) = 0;
    res.at<double>(1, 0) = sin(theta);
    res.at<double>(1, 1) = cos(theta);
    res.at<double>(1, 2) = 0;
    res.at<double>(2, 0) = 0;
    res.at<double>(2, 1) = 0;
    res.at<double>(2, 2) = 1;
    return res;
}

Point3f vector_rotation(float theta, Vec3f vect) {
    // manual rotation of vect by angle (in radians) theta along z axis
    auto co = cos(theta);
    auto si = sin(theta);
    auto vec_x = vect[0];
    auto vec_y = vect[1];
    //cout << vec_x << " " << vec_y << endl;
    Point3f res(co * vec_x - si * vec_y, co * vec_y + si * vec_x,vect[2]);
    return res;
}

Point3f new_vector_angle(double theta, Vec3f vect) {
    // rotate the vector s.t the angle with horizontal line is theta
    auto vec_x = vect[0];
    auto vec_y = vect[1];
    auto alpha = get_angle(vect);  // alpha is the current angle with the horizontal line
    //cout << "Theta, alpha: " << theta*rtd << "; " << alpha*rtd << endl;
    auto co = cos(-alpha+theta);      // rotate by -alpha
    auto si = sin(-alpha+theta);
    //cout << vec_x << " " << vec_y << endl;
    Point3f res(co * vec_x - si * vec_y, co * vec_y + si * vec_x, vect[2]);
    return res;
}

void set_vector_rotation(Skeleton3D* ske, Skeleton3D* child, double angle, string method = "rotation") {
    // choose "rotation" or "set angle" to either rotate the vector or set its position for a given angle
    auto grand_child_pos = child->get_root();
    Point3f temp;

    if (method == "rotation") {
        temp = vector_rotation(angle * dtr, grand_child_pos->getPos());
    }
    else if (method == "set angle"){
        temp = new_vector_angle(angle * dtr, grand_child_pos->getPos());
    }
    child->get_root()->setPos(temp);
    ske->updateAbsolutePosition();
}

/*Mat optimize(Mat& im, Skeleton2D& t){
    Mat c, c2, c3;
    c = t.toMat(im.rows, im.cols);
    
    ////Min Max Scaling using bboxes
    //Mat h_rect = find_rect_match(im, c);
    //t.transform(h_rect);
    //c2 = t.toMat(im.rows, im.cols);
    //
    ////Getting rot and trans using pca
    ////Mat h = get_rot(im, c2);
    //Mat h = get_trans(im, c2);
    //t.transform(h);
    //c3 = t.toMat(im.rows, im.cols);
    //
    return c;
}
*/

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

void show_2D_image(Skeleton3D skeleton, Mat image, bool show = false) {
    // set show = true to show the 2D skeleton (with names and arrows)
    Mat diag = Mat::eye(3, 3, CV_32FC1);
    skeleton.updateAbsolutePosition();
    Skeleton2D t = skeleton.project();
    Mat c = t.toMat(image.rows, image.cols, show);
    auto r = local_iou(image, c);
    show_merged("recuit", image, c, to_string(r[0])); waitKey(0);
}

double get_score(Skeleton3D ske, Mat image, bool show = false, string shape = "ellipse") {
    Skeleton2D t = ske.project_individual();
    if (ske.get_name() == "hip") {      // get the entire skeleton
        //ske.updateAbsolutePosition();
        t = ske.project();
        //cout << ske.get_child(0)->get_child(1)->get_child(0)->get_root()->getAbsPos() << endl;
    }
    //cout << ske.get_name() << " abs pos: " << ske.get_root()->getAbsPos() << endl;
    Mat c = t.toMat(image.rows, image.cols, show, shape);
    auto r = local_iou(image, c);
    if (show) show_merged("recuit", image, c, to_string(r[0])); waitKey(10);
    return r[0];
}



/*double try_and_optimize(Mat& target, Vec3f& v, Skeleton3D& s, double angle_min, double angle_max){ 
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
            imshow("c", c); waitKey(5);
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

        cout << "name, abs pos and pos: " << moving->get_name() << ", " << moving->get_root()->getQuat() << ", " << moving->get_root()->getPos() << endl;
        //cout << "v and angle: " << v << ", " << angle << endl;
        Quaternion q(angle, v);
        moving->rotate(q);
        cout << "name, abs pos and pos: " << moving->get_name() << ", " << moving->get_root()->getAbsPos() << ", " << moving->get_root()->getPos() << endl << endl;

        s.updateAbsolutePosition();
        
        Skeleton2D t = s.project(rep);
        //t.normalize(target.rows, target.cols);
        cout << t.get_name();
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
    //Mat rep = Mat::eye(3, 3, CV_32FC1);
    
    char key = 'a';
    ofstream fichier("energy.txt", ios::out | ios::trunc);  // creation of the file energy.txt
    
    if(fichier)
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
                cout << "r: " << r << endl;
                if (j==0){
                    k_b = r/(100*T);
                }
                //double alea = (double) rand()/RAND_MAX;
                if (r>r_max) {
                    fichier<<"0;";
                    r_max = r;
                    s_max = s;
                }
                else {
                    fichier<<"2;";
                    new_confs++;
                }
                fichier<<r_max<<";"<<T<<endl;
            }
            if (j%5==0){
                //cout<<"T="<<T<<endl;

                cout << "project" << rep << endl;
                Skeleton2D t = s_max.project(rep);
                //t.normalize(target.rows, target.cols);
                Mat c = optimize(target, t);
                show_merged("recuit", target, c, to_string(r_max));
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
*/

Skeleton3D* get_child(string name, map<string, vector<int>> m, Skeleton3D* ske) {
    // select child knowing its relative position to the root
    auto path = m[name];
    //Skeleton3D* res = ske.get_child(path[0]);
    //for (size_t i = 1; i < path.size();i++) {
    //    res = res->get_child(path[i]);
    //}
    if (path.empty()) {
        return ske;
    }
    Skeleton3D* res = ske->get_child(path[0]);
    path.erase(path.begin());
    m[name] = path;
    //return res;
    return get_child(name, m, res);
}

map<string, vector<int>> create_map(Skeleton3D* ske) {
    map<string, vector<int>> m;
    for (int i = 0; i < ske->get_children_size(); i++) {
        m[ske->get_child(i)->get_name()] = ske->get_child(i)->get_hierarchy();
        auto temp_m = create_map(ske->get_child(i));
        m.insert(temp_m.begin(), temp_m.end());
    }
    return m;
}

int get_direction(Skeleton3D* ske, Skeleton3D* bone, Mat image) {
    //cout << "angle: " << get_angle(limb->get_root()->getPos()) << endl;
    //cout << ske->get_child(0)->get_child(1)->get_child(0)->get_root()->getPos() << endl;

    //cout << "debug 0.001: " << get_score(*bone, image, false, "circle") << endl;
    if (get_score(*bone, image, false, "circle") == 1) {        // if it is already inside, go either way
        return 1;
    }

    float angle = 5;
    auto r = get_score(*bone, image, false);
    //get_score(*ske, image, true);

    set_vector_rotation(ske, bone, angle);
    auto r1 = get_score(*bone, image, false);
    //get_score(*ske, image, true);
    
    set_vector_rotation(ske, bone, -2*angle);
    auto r2 = get_score(*bone, image, false);
    //get_score(*ske, image, true);

    //cout << r1 << ", " << r << ", " << r2 << endl;

    set_vector_rotation(ske, bone, angle);          // go back to original position
    //return r1 < r2 ? "r" : "l";     // if r2 is a biggest score, should rotate counter clock-wise (i.e to the right as all vectors face down)
    return r1 < r2 ? -1 : 1;          // so should mult by -1
}

double get_border(Skeleton3D* ske, Skeleton3D* bone, Mat im, int signe, double angle){
    // goal is to find the two angles s.t the extremity of the bone is on either sides of the leg (or any other limb)
    // => take the average to get the middle of the leg
    if (angle < 1) {
        return get_angle(bone->get_root()->getPos()) * rtd;
    }
    //int signe = direction == "r" ? -1 : 1;              // if direction == r, angle is negative
    auto is_inside = get_score(*bone, im, false, "circle") < 1 ? 0 : 1;      // if is_inside == 1 => extremity is inside the shape : else 0
    auto min_angles = bone->get_min_angle_constraints();
    auto max_angles = bone->get_max_angle_constraints();
    //cout << "max angle bis: " << bone->get_max_angle_constraints() << ", " << bone->get_min_angle_constraints() << endl;
    auto temp = get_angle(bone->get_root()->getPos()) * double(rtd);
    while (get_score(*bone, im, true, "circle") == is_inside && min_angles <= temp && temp <= max_angles) {
        set_vector_rotation(ske, bone, signe * angle);
        temp = get_angle(bone->get_root()->getPos()) * double(rtd);
        //cout << "DEBUG: " << min_angles << ", " << temp << ", " << max_angles << endl;
    }
    if (get_score(*bone, im, false, "circle") == is_inside) {           // if it is still inside but reached to min/max angles
        return get_border(ske, bone, im, -signe, angle);                // go the other direction
    }
    return get_border(ske, bone, im, -signe, angle / 4);
}

void place_middle(Skeleton3D* ske, Skeleton3D* bone, Mat im) {
    // place the extremity of the bone in the center of the leg (or other limb)
    auto d = get_direction(ske, bone, im);
    auto r = get_score(*bone, im, false, "circle");

    //cout << "direction and score: " << d << ", " << r << endl;

    auto first_angle = get_border(ske, bone, im, d, 5);
    // auto n = log(10/2)/log(4);
    if (r == 1) d *= -1;
    set_vector_rotation(ske, bone, double(d) * 5);
    auto second_angle = get_border(ske, bone, im, d, 5);
    auto average = (second_angle + first_angle) / 2;
    /*cout << "name: " << bone->get_name() << endl;
    cout << "first angle: " << first_angle << endl;
    cout << "second angle: " << second_angle << endl;
    cout << "average: " << average << endl;
    cout << "max angle: " << bone->get_max_angle_constraints() << ", " << bone->get_min_angle_constraints() << endl << endl;
    */set_vector_rotation(ske, bone, average, "set angle");
    //show_2D_image(*ske, im, true);
    for (int i = 0; i < bone->get_children_size();i++) {
        place_middle(ske, bone->get_child(i), im);
    }
}

void optimisation_middle(Skeleton3D* ske, map<string, vector<int>> m, Mat im) {
    auto child = get_child("xrf", m, ske);
    place_middle(ske, child, im);

    child = get_child("xlf", m, ske);
    place_middle(ske, child, im);

    child = get_child("xlb", m, ske);
    place_middle(ske, child, im);

    child = get_child("xrb", m, ske);
    place_middle(ske, child, im);

    child = get_child("nose", m, ske);
    place_middle(ske, child, im);
}

vector<double> get_intersection_circle_line(double x0,double y0, double r, Skeleton3D lower_bone, Skeleton3D upper_bone) {
    // find x,y s.t 
    // (x-x0)² + (y_y0)² = r² and
    // y = a*x + b
    auto coord_lower = lower_bone.get_root()->getAbsPos();
    auto coord_upper = upper_bone.get_root()->getAbsPos();
    /*coord_lower.y = -6.64;
    coord_lower.x = -0.27;
    coord_upper.y = -5;
    coord_upper.x = 0;*/
    auto a = (coord_lower.y - coord_upper.y) / (coord_lower.x - coord_upper.x);
    auto b = coord_lower.y - a * coord_lower.x;
    auto root = sqrt((a * a + 1) * r * r + b * (2 * y0 - 2 * a * x0) - pow(y0 - a * x0, 2) - b * b);
    auto xp = (-a * b + a * y0 + x0 + root) / (a * a + 1);
    auto xm = (-a * b + a * y0 + x0 - root) / (a * a + 1);
    auto yp = a * xp + b;
    auto ym = a * xm + b;
    auto dist_p = sqrt(pow(xp - coord_upper.x, 2) + pow(yp - coord_upper.y, 2));
    auto dist_m = sqrt(pow(xm - coord_upper.x, 2) + pow(ym - coord_upper.y, 2));
    vector<double> p = { xp,yp, coord_upper.z};
    vector<double> m = { xm,ym, coord_upper.z};
    return dist_p < dist_m ? p : m;     // return the closest one
}

void update_child(Skeleton3D* ske, Skeleton3D* bone, Skeleton3D* child) {
    // bone is the new position of the upper part of the leg (i.e "rf", "lf" etc)
    // child is the child of bone before rotation
    auto bone_pos = bone->get_root()->getAbsPos();
    auto grand_child = child->get_child(0);
    auto child_pos = child->get_root()->getPos();
    double r = sqrt(pow(child_pos.x, 2) + pow(child_pos.y, 2));         // length of child
    auto temp = get_intersection_circle_line(bone_pos.x, bone_pos.y, r, *grand_child, *child);
    Point3f child_new_abs_pos(temp[0], temp[1], temp[2]);
    auto child_new_pos = child_new_abs_pos - bone_pos;
    cout << "DEBUG: " << sqrt(pow(child_new_pos.x, 2) + pow(child_new_pos.y, 2)) << endl;
    auto angle = get_angle(child_new_pos) * double(rtd);            // convert in degree
    set_vector_rotation(ske, bone->get_child(0), angle, "set angle");
}

double improve_legs(Skeleton3D* ske, Skeleton3D* bone, Mat im) {
    auto r = get_score(*bone->get_child(0), im, false, "ellipse");
    double r_max = 0;
    double angle;
    auto previous_pos = bone->get_root()->getPos();
    for (int i = 0; i < 10; i++) {
        cout << "previous 1: " << previous_pos << endl;
        set_vector_rotation(ske, bone, pow(-1, i) * 3 * i);                 //0, -3, 6 ...
        r = get_score(*bone->get_child(0), im, true, "ellipse");
        if (r > r_max) {
            r_max = r;
            angle = get_angle(bone->get_child(0)->get_root()->getPos());
        }
    }
    return angle;
}

int main(int argc, const char * argv[]) {
 
    int n = 7;
    std::string path_temp = argv[0];
    String file_path = path_temp + "/../../../cpp/data/elephants2_mask/";

    //SnapshotGraph g(file_path, n, true);
    //vector<int> path = g.get_path();

    vector<int> path = { 3,2,1 };           // use the two lines above to compute the order 

    Vec3f v_i(1, 0, 5);
    v_i = v_i/norm(v_i);
    Sprite s(file_path, path, false, v_i);
    char key = 'a';
    int theta = 30;
    /*while (key!='q'){
        s.display_one_horiz((double)theta);
        s.display_one_vert((double)theta); key = waitKey();
        if ((int)key==0)
            theta = (theta+5)%360;
        if ((int)key==1)
            theta = (theta-5)%360;
    }*/
    
    int i_query = 1;
    stringstream query;
    query << file_path << i_query << ".png";
    Mat im = imread(query.str(), cv::IMREAD_GRAYSCALE);
    
    if (i_query == 1) resize(im, im, Size(), 1.6, 1.6);
    if (i_query == 2) resize(im, im, Size(), 1.4, 1.4);
    if (i_query == 3) resize(im, im, Size(), 2, 2);
    //imshow("hello", im);waitKey(0);
    cout<<"Image: " << query.str() << ", with shape = "<<im.size<<endl;
    
    Skeleton3D hip = test3D();
    hip.updateAbsolutePosition();
    Vec3f v(0, 0, 1);
    v = v/norm(v); 
    Mat rep = repere_plane(v);

    Mat h = Mat::eye(3, 3, CV_32FC1);
    Mat diag = Mat::eye(3, 3, CV_32FC1);
    h.at<float>(0, 0) = 35;
    h.at<float>(1, 1) = 40;
    hip.transform(h);
    if (i_query == 1) hip.transform_translate(90, 90, 0);
    if (i_query == 2) hip.transform_translate(80, 75, 0);
    if (i_query == 3) hip.transform_translate(105, 100, 0);

    hip.create_hierarchy();
    hip.add_angle_constraints();
    map<string, vector<int>> m;
    m = create_map(&hip);       // map name of bone with its relative position to the root

    //auto child = get_child("lf", m, &hip);
    //cout << "name: " << child->get_name() << child->get_root()->getPos() << endl;

    /*for (auto i : m) {
        cout << "m: " << i.first;
        for (int j = 0; j < i.second.size(); j++) {
             cout << "; " << i.second[j];
        }
        cout << endl;
    }*/
    
    //auto grand_child = hip.get_child(0)->get_child(0);
    //cout << "hip grand child: " << grand_child->get_name() << grand_child->get_root()->getPos() << endl;

    show_2D_image(hip, im, true);

    optimisation_middle(&hip, m, im);

    show_2D_image(hip, im, true);

    Skeleton3D child = *get_child("xrf", m, &hip);
    cout << "child: " << child.get_root()->getAbsPos() << endl;

    auto bone = get_child("rf", m, &hip);

    set_vector_rotation(&hip, bone, -10);

    cout << child.get_root()->getAbsPos() << endl;

    update_child(&hip, bone, &child);
    //improve_legs(&hip, get_child("rf", m, &hip), im);

    show_2D_image(hip, im, true);

    /*auto lower_child = get_child("xxrf", m, &hip);
    auto upper_child = get_child("xrf", m, &hip);

    auto test = get_intersection_circle_line(1.84, -2.572, 2.236, *lower_child, *upper_child);*/

    //recuit_simule(im, v, hip);

    cout << "Key: " << key << endl;

    im.release();
    return 0;
}
