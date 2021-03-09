#include <iostream>
#include <iomanip>
#include "electronicmap/common.h"
#include "electronicmap/config.h"
#include "electronicmap/coordinates.h"

using namespace std;

int main(int argv, char** argc) {
    char* path = "./data/points.txt";
    DMS central(113, 0, 0);
    CConfigCoors config(path, central);
    config.SetCentral(central);
    CCoors coor(config);
    auto data = coor.BLH2XYZ_Batch();
    auto data_blh = coor.XYZ2BLH_Batch();
    auto blh_origin = coor.GetPointsBLH();

    for (auto da : data)
        cout << da.X_ << "  " << da.Y_ << "  " << da.Z_ << endl << endl;;

    for(int i = 0; i < data_blh.size(); ++i) 
        cout << setprecision(15) << data_blh[i].B_Rad_ - blh_origin[i].B_Rad_ << "   " <<
                data_blh[i].L_Rad_ - blh_origin[i].L_Rad_ << "   " <<
                data_blh[i].H_ - blh_origin[i].H_ << "   " << endl;
}