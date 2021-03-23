/*---------------------------------------------
   dijkstra.cpp
   create on 09 Mar 2021 ZHUOXU WHU
---------------------------------------------*/
#include "electronicmap/dijkstra.h"

CDijkstra::CDijkstra() {
    dist_ = NULL; isFind_ = NULL; path_ = NULL;
    mindist_ = 0;
}

CDijkstra::CDijkstra(CDijksConfig config) {
    adj_matrix_ = config.GetAdjMat();
    vertex_num_ = config.GetVerNum();
    isFind_ = new bool[vertex_num_];
    dist_ = new double[vertex_num_];
    path_ = new int[vertex_num_];
    memset(isFind_, 0, sizeof(bool) * vertex_num_);
    // memset(dist_, __DBL_MAX__, sizeof(double) * vertex_num_); // will be nan
    // will succeed
    for(int i = 0; i < vertex_num_; ++i) 
        dist_[i] = __DBL_MAX__;
    memset(path_, -1, sizeof(int) * vertex_num_);
    mindist_ = 0;
}

CDijkstra::~CDijkstra() {
    if(isFind_ != NULL && dist_ != NULL && path_ != NULL) {
        delete[] isFind_;
        delete[] dist_;
        delete[] path_;
    }
    dist_ = NULL; isFind_ = NULL; path_ = NULL;
}

VecInt CDijkstra::FindPathNear(int verNumStart, int verNumEnd) {
    // legal check
    VecInt result;
    if(verNumStart >= vertex_num_ || verNumEnd >= vertex_num_) {
        cout << "please check your input, vertex index wrong" << endl;
        return result;
    }
    GenerateMap(verNumStart);
    for(int i = verNumEnd; i != verNumStart; i = path_[i])
        result.push_back(i);
    result.push_back(verNumStart);
    std::reverse(result.begin(), result.end());
    ofstream out("./data/path.txt", ios::app);
    out << "from " << verNumStart << " to " << verNumEnd << " path is " << endl;
    for(auto re : result) out << re << " ";
    out << endl << "distance is " << dist_[verNumEnd] << endl << endl;
    out.close();
    return result;
}

int CDijkstra::GenerateMap(int verNumStart) {
    int next = verNumStart;
    int prev = next;
    isFind_[next] = 1; dist_[next] = 0; 
    ofstream out("./data/result.txt", ios::app);
    out << "起点为 " << next << endl;
    while (true) {
        out << "当前dist数组为： " << endl;
        for(int i = 0; i < vertex_num_; ++ i) out << dist_[i] << " ";
        out << endl;
        out << "当前path数组为： " << endl;
        for(int i = 0; i < vertex_num_; ++ i) out << path_[i] << " ";
        out << endl;
        out << "当前S数组为： " << endl;
        for(int i = 0; i < vertex_num_; ++ i) out << isFind_[i] << " ";
        out << endl;
        auto attach = SearchAttach(next, prev);
        UpdateState(next, attach, mindist_);
        out <<"查询到连接于 " << next << "的顶点有： " << endl;
        for(auto da : attach) out << da << "  ";
        out << "距离分别为： " << endl;
        for(auto da : attach) out << adj_matrix_(next, da) << "  ";
        out << endl;
        prev = next;
        next = DistMin();
        if (next == -1) 
            break;
        out << "当前距离 " << prev << " 最近的点为 " << next << " 距离是" << mindist_ << endl;
        isFind_[next] = true;
        // cout << next << endl;
        // for(int i = 0; i < vertex_num_; ++i)
        //     cout << dist_[i] << "   " << path_[i] << "  " << isFind_[i] << endl;
        out << endl;
    }
    out << endl;
    out.close();
}


VecInt CDijkstra::SearchAttach(int vernum, int startNum) {
    VecInt result;
    for (int j = 0; j < vertex_num_; ++ j) {
        if(vernum == j || adj_matrix_(vernum, j) == __DBL_MAX__)
            continue;
        if(j == startNum) continue;
        result.push_back(j);
    }
    return result;
}

void CDijkstra::UpdateState(int vertex, VecInt attaches, double mindist) {
    for(auto attach : attaches) {
        if(isFind_[attach] == false) {
            dist_[attach] = mindist + adj_matrix_(vertex, attach);
            path_[attach] = vertex;
        }
    }
}

int CDijkstra::DistMin() {
    double min = __DBL_MAX__;
    int index = -1;
    for(int i = 0; i < vertex_num_; ++i)
        if(dist_[i] < min && isFind_[i] == false) {
            index = i; min = dist_[i];
        }
    mindist_ = min;
    return index;
}