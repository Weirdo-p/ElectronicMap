/*---------------------------------------------
   coorconfig.cpp
   create on 09 Mar 2021 ZHUOXU WHU
---------------------------------------------*/
#include "electronicmap/coorconfig.h"
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

CConfigCoors::CConfigCoors() { }

CConfigCoors::CConfigCoors(char* FilePath, DMS central) {
    double central_deg = 0;
    central_dms_ = central;
    elli_ = WGS84;

    DMS2Deg(central_dms_, central_deg);
    Deg2Rad(central_deg, central_rad_);
    ReadFile(FilePath);
}

bool CConfigCoors::ReadFile(char* FilePath) {
    ifstream in(FilePath);
    if(!in) {
        cerr << "can not open file" << endl;
        exit(1);
    }

    while (!in.eof()) {
        string line; 
        stringstream ss;
        DMS B, L;
        double H;
        int num;
        getline(in, line);
        if(line[0] == '#')
            continue;
        ss << line;
        ss >> num >> B.Deg_ >> B.Min_ >> B.Sec_ >> L.Deg_ >> L.Min_ >> L.Sec_ >> H;
        points_.push_back(BLH(B, L, H));
    }
    elli_= WGS84;

    return true;
}

void CConfigCoors::SetCentral(DMS central) {
    double deg = 0;
    central_dms_ = central_dms_;
    DMS2Deg(central_dms_, deg);
    Deg2Rad(deg, central_rad_);
}

void CConfigCoors::SetEllipsoid(EllipsoidType type) {
    elli_ = type;
}

DMS CConfigCoors::GetCentralDms() {
    return central_dms_;
}

double CConfigCoors::GetCentralRad() {
    return central_rad_;
}

vector<BLH> CConfigCoors::GetPoints() {
    return points_;
}

EllipsoidType CConfigCoors::GetEllipsoidType() {
    return elli_;
}
