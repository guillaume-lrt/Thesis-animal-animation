//
//  animating.cpp
//  projetc
//
//  Created by Coralie DAVID on 17/12/2017.
//  Copyright © 2017 Manon Romain. All rights reserved.
//

#include "animating.hpp"


SnapshotGraph::SnapshotGraph(String path, int n, bool cycle){
    String file = path+"graph.txt";
    try{
        load_from_txt(file);
        
    }
    catch (Exception e){
        e.formatMessage();
        full_cycle = cycle;
        ShapeContextDistance scd;
        Size s = SIZE_I;
        vector<Mat> vm;
        for (int i_query = 0; i_query<n; i_query++){
            vector<double> gi(i_query);
            stringstream query;
            query << path << i_query << ".png";
            Mat im = imread(query.str(), cv::IMREAD_GRAYSCALE);
            cout<< "Image size" <<im.size()<<endl;
            resize(im, im, s);
            vm.push_back(im);
        }
        g = scd.apply(vm);
        write_in_file(file);
    }
    delete_outliers(2.0);
    select_min_max();
    auto t = chrono::system_clock::now();
    create_path();
    auto t2 = chrono::system_clock::now();
    chrono::duration<double> d = t2 - t;
    cout<<"Time took "<<d.count()<<"s"<<endl;
}

void SnapshotGraph::load_from_txt(String file){

    ifstream fichier(file, ios::in);
    
    if(fichier)  
    {
        string ligne;
        getline(fichier, ligne);
        int n = stoi(ligne);
        g = vector<vector<double>>(n);
        for (int i = 0; i<n; i++)
            g[i] = vector<double>(n);
        while(getline(fichier, ligne))
        {
            vector<String> v = explode(ligne, ';');
            g[stoi(v[0])][stoi(v[1])] = stod(v[2]);
        }

        fichier.close();
    }
    else {
        cerr << "Fichier text non présent !" << endl;
        throw Exception();
    }
    
}

void SnapshotGraph::write_in_file(String file){
    
    ofstream fichier(file, ios::out | ios::trunc);
    
    if(fichier)
    {
        int n = g.size();
        fichier<<n<<endl;
        for (int i = 0; i<n; i++){
            for (int j = 0; j<n; j++){
                fichier<<i<<";"<<j<<";"<<g[i][j]<<endl;
            }
        }
        fichier.close();
    }
    else
        cerr << "Erreur à l'ouverture pour écrire !" << endl;
}


void SnapshotGraph::delete_outliers(double c){
    /** remove the nodes if the distance with the closest node is greater than some Gaussian threshold
    i.e if min_dist > µ+c∗σ with c = 2.0 **/
    size_t n = g.size();
    vector<double> min_dist(n);
    
    //Compute min dist to graph
    min_dist[0] = g[1][0];
    for (int i = 0; i<n; i++){
        if (i>0)
            min_dist[i] = g[i][0];
        vector<int> idx(n);
        
        for (int ii=0; ii<n; ++ii){
            if (i==ii) continue;
            double d = g[i][ii];
            min_dist[i] = min(d, min_dist[i]);
        }
        
    }
    
    vector<bool> b(n);
    double mu = avg(g);
    cout<<mu<<endl;
    double std = sqrt(variance(g, mu));
    cout<<"std:"<<std<<endl;
    for (int i = 0; i<n; i++){
        b[i] = (min_dist[i]-mu>c*std);
        if (verbose){
            cout<<"Min dist of "<<i<<" is "<<min_dist[i]<<endl;
            if (b[i])
                cout<<i<<" is an outlier"<<endl;
            else
                cout<<i<<" isn't an outlier"<<endl;
        }
    }
    this->outliers = b;
}

void SnapshotGraph::select_min_max(){
    double dmax = g[1][0];
    start = 1; end = 0;
    for (int i = 1; i<g.size(); i++){
        for (int j=0; j<i; j++){
            if (g[i][j]>alpha){
                alpha = g[i][j];            //alpha = max_{k,l}(D(S_k,S_l))  i.e max distance between two nodes in the graph
            }
            if (outliers[j] || outliers[i])
                continue;
            if (g[i][j]>dmax){
                dmax = g[i][j];
                start = i;
                end = j;
            }
        }
    }
    cout<<"start is "<<start<<endl;
    cout<<"end is "<<end<<endl;
};



vector<int> SnapshotGraph::generate_path_half_cycle(int l){
    vector<int> myvector;
    for (int i=0; i<g.size(); ++i){
        if (i==start || i==end) continue;
        myvector.push_back(i);
    }
    random_shuffle(myvector.begin(), myvector.end());
    vector<int> len_l;
    len_l.push_back(start);
    for (int i=0; i<l-2; ++i) len_l.push_back(myvector[i]);
    len_l.push_back(end);
    return len_l;
}

vector<int> SnapshotGraph::generate_path_full_cycle(int l){
    vector<int> myvector;
    for (int i=0; i<g.size(); ++i){
        if (i==start || i==end) continue;
        myvector.push_back(i);
    }
    random_shuffle(myvector.begin(), myvector.end());
    
    vector<int> len_l;
    len_l.push_back(start);
    int pos_end = rand()%(l-2);
    for (int i=0; i<l-1; ++i){
        if (i==pos_end) len_l.push_back(end);
        else len_l.push_back(myvector[i]);
    }
    len_l.push_back(start);
    return len_l;
}

vector<int> SnapshotGraph::generate_path(int l){
    if (full_cycle)
        return generate_path_full_cycle(l);
    else
        return generate_path_half_cycle(l);
}

double SnapshotGraph::energy_path(int l, vector<int> p){
    /** cost function that measure the difference between images
    Goal is to minimize this function **/
    vector<double> dist;
    double Cs = 0, Cu = 0, Cd = 0;
    // Cs measures local similarity => mean of distances between two consecutives nodes

    
    for (int i=0; i<l; ++i){
        if (i<l-1)
            dist.push_back(g[p[i]][p[i+1]]);
        for (int j = 0; j<l; ++j){
            if (j==i) continue;
            Cd += g[p[i]][p[j]];
        }
    }
    Cs = avg(dist);
    Cu = variance(dist, Cs);
    
    Cd *= alpha/(l*(l-1));
    Cs *= alpha;
    return Cs + lambda_1*Cu + lambda_2*(1-Cd);
}


void SnapshotGraph::create_path(){
    /***
     Simulated annealing - it don't understand why
     ***/
    double k = 500;
    double T = 50000;
    double limit = 0.01;
    
    int min_L;
    
    if (full_cycle)
        min_L = 3;
    else
        min_L = 2;
    
    int L = (rand()%(g.size()-min_L+1))+min_L;
    if (L==g.size()) cout<<"DEBUG"<<endl;
    vector<int> curr_path = generate_path(L);
    double E = energy_path(L, curr_path);
    
    int new_L;
    vector<int> new_path;
    double new_E, delta_E;
    
    while (T>limit){
        for (int i=0; i<k; ++i){
            new_L = rand()%(g.size()-min_L) + min_L;
            new_path =  generate_path(new_L);
            new_E = energy_path(new_L, new_path);
            delta_E = new_E - E;
            if ((delta_E<=0)){
                E = new_E;
                L = new_L;
                curr_path = new_path;
            }
        }
        T *= 0.992;
    }
    cout<<"TERMINE"<<endl;
    E = energy_path(L, curr_path);
    cout<<"E is "<<E<<endl;
    for (int i: curr_path)
        cout<<i<<";";
    cout<<endl;
    path = curr_path;
};



