#include <iostream>
#include <iomanip>
#include "electronicmap/common.h"
#include "electronicmap/config.h"

using namespace std;

int main(int argv, char** argc) {
    Ellipsoid wgs84(WGS84);
    cout << setprecision(15) << wgs84.a_ << endl <<
            wgs84.b_ << endl <<
            wgs84.c_ << endl << 
            wgs84.e2_ << endl <<
            wgs84.eprime2_ << endl <<
            wgs84.alpha_ << endl << endl;
    double deg = 180, rad = 0;
    Deg2Rad(180, rad);
    cout << rad << endl << endl;

    DMS B(360, 0, 0), L(114, 0, 0);
    BLH blh(B, L, 0);
    cout << blh.B_Rad_ << endl <<
            blh.L_Rad_ << endl << endl;
    char* path = "./data/points.txt";
    DMS central(113, 0, 0);
    CConfigCoors config(path, central);
    config.SetCentral(central);
}