/*---------------------------------------------
   dijksconfig.cpp
   create on 09 Mar 2021 ZHUOXU WHU
---------------------------------------------*/

#include "electronicmap/dijkscommon.h"
#include "electronicmap/dijksconfig.h"

#include <fstream>
#include <sstream>

CDijksConfig::CDijksConfig() { }

CDijksConfig::CDijksConfig(char* path) {
    ReadFile(path);
}

CDijksConfig::CDijksConfig(MatNNd adj_matrix) {
    adj_matrix_ = adj_matrix;
}

MatNNd CDijksConfig::GetAdjMat() {
    return adj_matrix_;
}


void CDijksConfig::ReadFile(char* filepath) {

    ifstream in(filepath);
    if(!in) {
        cerr << "can not open file" << endl;
        exit(-1);
    }
    while (!in.eof()) {
        string line;
        getline(in, line);
        if (line.length() == 1) {
            InitAdjMatrix(line);
            continue;
        }
        SetAdjMatrix(line);
    }
    for (int i = 0; i < adj_matrix_.rows(); ++i)
        for (int j = 0; j < adj_matrix_.cols(); ++j)
            if (i == j)
                adj_matrix_(i, j) = 0;
            else if (i != j && adj_matrix_(i, j) == 0)
                adj_matrix_(i, j) = __DBL_MAX__;
}

void CDijksConfig::InitAdjMatrix(string line) {
    stringstream buff;
    int num = 0;
    buff << line;
    buff >> num;
    if(num <= 0) {
        cerr << "vertex num error" << endl;
        exit(-1);
    }
    vertex_num_ = num;
    adj_matrix_.resize(num, num);
}

int CDijksConfig::GetVerNum() {
    return vertex_num_;
}


void CDijksConfig::SetAdjMatrix(string line) {
    stringstream buff;
    buff << line;
    int start, end, weight;
    buff >> start >> end >> weight;
    if (weight < 0) {
        cerr << "error! minus weight can not use dijstra algorithm" << endl;
        exit(-1);
    }
    try {
        adj_matrix_(start - 1, end - 1) = weight;
    }
    catch(...) {
        cerr << "error happened when set weights to adjacency matrix" << endl;
        exit(-1);
    }
}
