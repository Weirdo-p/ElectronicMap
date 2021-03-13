/*---------------------------------------------
   mapproj.cpp
   create on 09 Mar 2021 ZHUOXU WHU
---------------------------------------------*/

#include "electronicmap/mapproj.h"

GaussProj::GaussProj() {
    elli_.SetEllipsoidParam(WGS84);
}

GaussProj::GaussProj(CConfigCoors config) {
    elli_.SetEllipsoidParam(config.GetEllipsoidType());
    central_dms_ = config.GetCentralDms();
    central_rad_ = config.GetCentralRad();
    points_blh_ = config.GetPoints();
    coor.SetEllipsoid(config.GetEllipsoidType());
}

bool GaussProj::GetN(double B,double &N) {
    double sinB = sin(B);
    double sub = sqrt(1 - elli_.e2_ * sinB * sinB);
    if(sub <= 1e-6)
        return false;
    N = elli_.a_ / sub;

    return true;
}
int GaussProj::GetL0(DMS longi) {
    double deg = longi.Deg_ + longi.Min_ / 60.0 + longi.Sec_ / 3600.0;
    int n = int(deg / 6);
    int L0 = 6 * n - 3;
    return L0;
}

double GaussProj::Geteta(double B) {
    double cosB = cos(B);
    double cosB_sqr = cosB * cosB;
    double eta = (elli_.eprime2_ * cosB_sqr);
    return eta;
}

double* GaussProj::GetCoef() {
    double* coef, a = elli_.a_, e2 = elli_.e2_;
    double e4 = e2 * e2, e6 = e4 * e2, e8 = e4 * e4;
    double e10 = e8 * e2;
    coef = new double[5];
    coef[0] = 1.0 + 3.0 / 4.0 * e2 + 45.0 / 64.0 * e4 + 175.0 / 256.0 * e6 + 11025.0 / 16384.0 * e8 + 43659.0 / 65536.0 * e10;
    coef[1] =       3.0 / 4.0 * e2 + 15.0 / 16.0 * e4 + 525.0 / 512.0 * e6 + 2205.0 / 2048.0 * e8 + 72765.0 / 65536.0 * e10;
    coef[2] =                        15.0 / 64.0 * e4 + 105.0 / 256.0 * e6 + 2205.0 / 4096.0 * e8 + 10395.0 / 16384.0 * e10;
    coef[3] =                                            35.0 / 512.0 * e6 + 315.0 / 2048.0 * e8 + 31185.0 / 131072.0 * e10;
    coef[4] =                                                                315.0 / 16384.0 * e8 + 3645.0 / 65536.0 * e10;

    return coef;
}

double GaussProj::Getx(double B, double t, double eta,
                       double l, double N, double X) {
    double coef[4] = {0}, sinB = sin(B), cosB = cos(B);
    coef[0] = X;
    coef[1] = N / 2.0 * sinB * cosB * pow(l, 2);
    coef[2] = N / 24.0 * sinB * pow(cosB, 3) * (5 - pow(t, 2) + 9 * pow(eta, 2) +
              4 * pow(eta, 4)) * pow(l, 4);
    coef[3] = N / 720.0 * sinB * pow(cosB, 5) * (61 - 58 * pow(t, 2) + pow(t, 4)) * pow(l, 6);
    double x = 0;
    for(int i = 0; i < 4; ++i)
        x += coef[i];
    return x;
}
double GaussProj::Gety(double B, double t, double eta, double l, double N) {
    double cosB = cos(B), coef[3] = {0};
    coef[0] = N * cosB * l;
    coef[1] = N / 6.0 * pow(cosB, 3) * (1 - pow(t, 2) + pow(eta, 2)) * pow(l, 3);
    coef[2] = N / 120.0 * pow(cosB, 5) * (5 - 18 * pow(t ,2) + pow(t, 4) + 14 *
              pow(eta, 2) - 58 * pow(eta, 2) * pow(t, 2)) * pow(l, 5);
    double y = 0;
    for (int i = 0; i < 3; ++i)
        y += coef[i];

    return y;
}

vector<XYZ> GaussProj::GetPointsPlain() {
    return points_plain_;
}

double GaussProj::GetBf(double x) {
    double* a = GetCoef();
    double tmp, Bf = x / (a[0] * elli_.a_ * (1.0 - elli_.e2_));
    tmp = Bf;
    double thres = 1e-9;
    int count = 0;
    while(true) {
        double coe = 0;
        for(int i = 1; i < 5; ++i)
            coe += a[i] / (i * 2) * sin(Bf * i * 2) * pow(-1, i + 1);
        double tempB = coe / a[0];
        tempB += tmp;
        if(abs(tempB - Bf) < thres)
            break;
        Bf = tempB;
        count++;
    }
    delete[] a;
    return Bf;
}

double GaussProj::GetM(double B) {
    double a = elli_.a_;
    double e2 = elli_.e2_;
    double sinB = sin(B);
    double sin2B = pow(sinB, 2);
    double supper = a * (1 - e2);
    double sub = pow((1.0 - e2 * sin2B), 1.5);
    double M = supper / sub;

    return M;
}

double GaussProj::GetB(double Bf, double y, double tf,
                       double Mf, double Nf, double etaf) {
    double etaf_sqr = etaf * etaf, tf_sqr = tf * tf, y2 = y * y;
    double Nf2 = Nf * Nf, Nf4 = Nf2 * Nf2;
    double coe[4], tmp = -y2 / (2.0 * Mf * Nf) * tf;
    coe[0] = Bf;
    coe[1] = 1 * tmp;
    coe[2] = -(y2 / (12.0 * Nf2) * (5.0 + etaf_sqr + 3.0 * tf_sqr - 9.0 * etaf_sqr * tf_sqr)) * tmp;
    coe[3] = tmp * (y2 * y2 / (360.0 * Nf4) * (61.0 + 90.0 * tf_sqr + 45.0 * tf_sqr));
    double B = 0;
    for(int i = 0; i < 4; ++i)
        B += coe[i];
    return B;
}

double GaussProj::GetL(double y, double Bf, double Nf, double tf, double etaf) {
    double etaf_sqr = etaf * etaf, tf_sqr = tf * tf;
    double cosBf = cos(Bf), y2 = y * y, Nf2= Nf * Nf;
    double coe[3], temp = y / (Nf * cosBf);
    coe[0] = 1.0 * temp;
    coe[1] = -y2 / (6.0 * Nf2) * (1.0 + etaf_sqr + 2.0 * tf_sqr) * temp;
    coe[2] = y2 * y2 / (120.0 * Nf2 * Nf2) * (5.0 + 6.0 * etaf_sqr + 28.0 * tf_sqr + 8.0 * etaf_sqr * tf_sqr + 24.0 * tf_sqr * tf_sqr) * temp;
    double l = 0;
    for(int i = 0; i < 3; ++i)
        l +=coe[i];
    l += central_rad_;
    return l;
}

XYZ GaussProj::Proj(BLH blh) {
    double N = 0;
    double B = blh.B_Rad_;
    GetN(B, N);
    double t = tan(B);
    double L0_rad = central_rad_;
    double l = blh.L_Rad_ - L0_rad;
    double eta = Geteta(B);
    double* a = GetCoef();
    double X = a[0] * B;
    for (int i = 1; i < 5; ++i)
        X += a[i] * sin(i * 2.0 * B) / (i * 2.0) * pow(-1, i);
    X *= (elli_.a_ * (1 - elli_.e2_));
    double x = Getx(B, t, eta, l, N, X);
    double y = Gety(B, t, eta, l, N);
    delete[] a;
    return XYZ(x, y, 0);
}

bool GaussProj::Proj_Batch() {
    for(auto point : points_blh_)
        points_plain_.push_back(Proj(point));
    return true;
}

BLH GaussProj::BackProj(XYZ xyz) {
    double Bf = GetBf(xyz.X_);
    double Mf = GetM(Bf);
    double tf = tan(Bf);
    double Nf = 0;
    GetN(Bf, Nf);
    double etaf = Geteta(Bf);
    double B = GetB(Bf, xyz.Y_, tf, Mf, Nf, etaf);
    double L = GetL(xyz.Y_, Bf, Nf, tf, etaf);
    return BLH(B, L, 0);    
}

vector<BLH> GaussProj::BackProj_Batch() {
    vector<BLH> result;
    for (auto pt : points_plain_)
        result.push_back(BackProj(pt));
    points_blh_ = result;
    return result;
}
