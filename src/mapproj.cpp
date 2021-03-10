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
    double eta = sqrt(elli_.eprime2_ * cosB_sqr);
    return eta;
}

double* GaussProj::GetCoef() {
    double* coef, a = elli_.a_, e2 = elli_.e2_;
    coef = new double[5];
    double m0 = a * (1 - e2);
    double m2 = 3.0 / 2.0 * e2 * m0;
    double m4 = 5.0 / 4.0 * e2 * m2;
    double m6 = 7.0 / 6.0 * e2 * m4;
    double m8 = 9.0 / 8.0 * e2 * m6;

    double a8 = m8 / 128.0;
    double a6 = m6 / 32.0 + m8 / 16.0;
    double a4 = m4 / 8.0 + 3.0 / 16.0 * m6 + 7.0 / 32.0 * m8;
    double a2 = m2 / 2.0 + m4 / 2.0 + 15.0 / 32.0 * m6 + 7.0 / 16.0 * m8;
    double a0 = m0 + m2 / 2.0 + 3.0 / 8.0 * m4 + 5.0 / 16.0 * m6 + 
                35.0 / 128.0 * m8;
    coef[0] = a0; coef[1] = a2; coef[2] = a4; coef[3] = a6; coef[4] = a8;
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
    double Bf = x / a[0];
    double thres = 1e-9;
    int count = 0;
    while(true) {
        double coe = 0;
        for(int i = 1; i < 5; ++i)
            coe += a[i] / (i * 2) * sin(Bf * i * 2) * pow(-1, i + 1);
        coe += x;
        double tempB = coe / a[0];
        if(abs(tempB - Bf) < thres)
            break;
        Bf = tempB;
        count++;
    }
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
    double etaf_sqr = etaf * etaf, tf_sqr = tf * tf;
    double coe[4];
    coe[0] = Bf;
    coe[1] = -tf / (2 * Mf * Nf) * y * y;
    coe[2] = tf / (24 * Mf * pow(Nf, 3)) * (5 + 3 * pow(tf, 3) + etaf_sqr - 9 * etaf_sqr * tf_sqr);
    coe[3] = -tf / (720 * Mf * pow(Nf, 5)) * (61 + 90 * tf_sqr + 45 * tf_sqr * tf_sqr) * pow(y, 6);
    double B = 0;
    for(int i = 0; i < 4; ++i)
        B += coe[i];
    return B;
}

double GaussProj::GetL(double y, double Bf, double Nf, double tf, double etaf) {
    double etaf_sqr = etaf * etaf, tf_sqr = tf * tf;
    double cosBf = cos(Bf);
    double coe[3];
    coe[0] = 1.0 / (Nf * cosBf) * y;
    coe[1] = -1.0 / (6 * pow(Nf, 3) * cosBf) * (1 + 2 * tf_sqr + etaf_sqr) * pow(y, 3);
    coe[2] = 1.0 / (120 * pow(Nf, 5) * cosBf) * (5 + 28 * tf_sqr + 24 * tf_sqr *
             tf_sqr + 6 * etaf_sqr + 8 * etaf_sqr * tf_sqr) * pow(y, 5);
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
