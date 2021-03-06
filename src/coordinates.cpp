/*---------------------------------------------
   coordinates.cpp
   create on 09 Mar 2021 ZHUOXU WHU
---------------------------------------------*/

#include "electronicmap/coordinates.h"

CCoors::CCoors() { }

CCoors::CCoors(CConfigCoors config) {
    points_blh_ = config.GetPoints();
    elli_.SetEllipsoidParam(config.GetEllipsoidType());
}

vector<XYZ> CCoors::BLH2XYZ_Batch() {
    XYZ xyz;
    vector<XYZ> result;
    for (auto point : points_blh_) {
        if(!BLH2XYZ(point, xyz)) {
            cout << "warning: error occurred when transforming BLH -> XYZ, " <<
                     "please check your input" << endl;
            break;
        }
        result.push_back(xyz);
    }
    points_xyz_ = result;
    return result;
}

bool CCoors::BLH2XYZ(BLH blh, XYZ& xyz) {
    if((!CheckB(blh.B_Dms_)) || (!CheckL(blh.L_Dms_)))
        return false;

    double sinB = sin(blh.B_Rad_);
    double cosB = cos(blh.B_Rad_);
    double sinL = sin(blh.L_Rad_);
    double cosL = cos(blh.L_Rad_);

    double N = 0;
    if(!GetN(blh.B_Rad_, N))
        return false;
    
    xyz.X_ = (N + blh.H_) * cosB * cosL;
    xyz.Y_ = (N + blh.H_) * cosB * sinL;
    xyz.Z_ = (N * (1 - elli_.e2_) + blh.H_) * sinB;

    return true;
}

bool CCoors::XYZ2BLH(XYZ xyz, BLH &blh) {
    double L_rad = atan2(xyz.Y_, xyz.X_);
    // 迭代初值  度为单位
    unsigned short int iteration = 0;
    // in deg
    double B0_deg = 1, B0_deg_iter = 1;
    double B0_rad;
    while(iteration != 15) {
        Deg2Rad(B0_deg, B0_rad);
        double sinB = sin(B0_rad);
        double N = 0;
        if(!GetN(B0_rad, N)) {
            cout << "error occurred when transforming XYZ -> BLH" << endl;
            break;
        }
        double H = xyz.Z_ / sin(B0_rad) - N * (1 - elli_.e2_);
        double up = xyz.Z_ + N * elli_.e2_ * sinB;
        double down = sqrt(xyz.X_ * xyz.X_ + xyz.Y_ * xyz.Y_);
        B0_rad = atan(up / down);
        Rad2Deg(B0_rad, B0_deg);
        blh.H_ = H;
        double error = abs(B0_deg - B0_deg_iter);
        if(error < 1e-20)
            break;
        B0_deg_iter = B0_deg;
        iteration ++;
    }
    blh.L_Rad_ = L_rad; blh.B_Rad_ = B0_rad;
    Rad2DMS(blh.L_Rad_, blh.L_Dms_); Rad2DMS(blh.B_Rad_, blh.B_Dms_);
    return true;
}
vector<BLH> CCoors::XYZ2BLH_Batch() {
    BLH blh;
    vector<BLH> result;
    for(auto point : points_xyz_) {
        XYZ2BLH(point, blh);
        result.push_back(blh);
    }
    points_blh_ = result;
    return result;
} 

bool CCoors::GetN(double B,double &N) {
    double sinB = sin(B);
    double sub = sqrt(1 - elli_.e2_ * sinB * sinB);
    if(sub <= 1e-6)
        return false;
    N = elli_.a_ / sub;

    return true;
}

void CCoors::SetEllipsoid(EllipsoidType type) {
    elli_.SetEllipsoidParam(type);
}

void CCoors::SetPointsBLH(vector<BLH> points_blh) {
    points_blh_ = points_blh;
}

void CCoors::SetPointsXYZ(vector<XYZ> points_xyz) {
    points_xyz_ = points_xyz;
}

vector<BLH> CCoors::GetPointsBLH() {
    return points_blh_;
}

vector<XYZ> CCoors::GetPointsXYZ() {
    return points_xyz_;
}

