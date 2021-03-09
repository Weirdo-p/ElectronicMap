#include "electronicmap/common.h"

DMS::DMS() {
    Deg_ = Min_ = Sec_ = 0;
}

DMS::DMS(int deg, int min, double sec) {
    Deg_ = deg; Min_ = min; Sec_ = sec;
}

Giant::Giant() {
    Integer_ = Frac_ = 0;
}

Giant::Giant(int Inte, double Frac) {
    Integer_ = Inte; Frac_ = Frac;
}

BLH::BLH() {
    B_Rad_ = L_Rad_ = H_ = 0;
}

XYZ::XYZ() {
    X_ = Y_ = Z_ = 0;
}

XYZ::XYZ(double x, double y, double z) {
    X_ = x; Y_ = y; Z_ = z;
}

BLH::BLH(DMS B, DMS L, double H) {
    double B_deg, L_deg;
    B_Dms_ = B; L_Dms_ = L; H_ = H;
    DMS2Deg(B_Dms_, B_deg); DMS2Deg(L_Dms_, L_deg);
    Deg2Rad(B_deg, B_Rad_); Deg2Rad(L_deg, L_Rad_);
}

BLH::BLH(double B, double L, double H) {
    B_Rad_ = B; L_Rad_ = L; H_ = H;
}

Ellipsoid::Ellipsoid() {
    type_ = WGS84;
    SetEllipsoidParam(type_);
}

Ellipsoid::Ellipsoid(EllipsoidType type) {
    type_ = type;
    SetEllipsoidParam(type);
}

bool Ellipsoid::SetEllipsoidParam(EllipsoidType type) {
    switch (type)
    {
    case WGS84: {
        a_ = 6378137.0;
        b_ = 6356752.3142;
        CalculateParam(a_, b_);
        break;
    }
    case CGCS2000: {
        a_ = 6378137.0;
        b_ = 6356752.3141;
        CalculateParam(a_, b_);
        break;
    }
    default: break;
    }
}

bool Ellipsoid::CalculateParam(double a, double b) {
    if (a == 0 || b == 0) return false;
    double e, e_prime;
    c_ = a * a / b;
    alpha_ = (a - b) / a;
    e = sqrt(a * a - b * b) / a;
    e_prime = sqrt(a * a - b * b) / b;
    e2_ = e * e;
    eprime2_ = e_prime * e_prime;

    return true;
}

void DMS2Deg(DMS dms, double &deg) {
    deg = 0;    // initialize
    deg += dms.Sec_ / 3600.0 + dms.Min_ / 60.0;
    if(dms.Deg_ < 0) {
        deg += (-dms.Deg_);
        deg = -deg;
    } else  deg += dms.Deg_;
}

void Deg2Rad(double deg, double &rad) {
    rad = (deg * PI) / 180.0;
}