/*---------------------------------------------
   coormain.cpp
   create on 09 Mar 2021 ZHUOXU WHU
---------------------------------------------*/
#include <iostream>
#include <iomanip>
#include "electronicmap/coorcommon.h"
#include "electronicmap/coorconfig.h"
#include "electronicmap/coordinates.h"
#include "electronicmap/mapproj.h"

using namespace std;

int main(int argv, char** argc) {
	// BLH->XYZ->BLH test
	if(argv != 3) {
		cout << "usage: coor [filepath] [certral longitude]" << endl;
		return 0;
	}
    char* path = argc[1];
    DMS central(atoi(argc[2]), 0, 0);
    CConfigCoors config(path, central);
    config.SetCentral(central);
    CCoors coor(config);
    auto data = coor.BLH2XYZ_Batch();
	for (auto da : data) cout << fixed << setprecision(4) << da.X_ << "  " << da.Y_ << "  " << da.Z_ << endl;
	cout << endl;
	coor.SetPointsXYZ(data);
	auto blhs_new = coor.XYZ2BLH_Batch();
	auto blhs = config.GetPoints();
	for(int i = 0; i < blhs.size(); ++i) {
		double diff_B = blhs_new[i].B_Rad_ - blhs[i].B_Rad_;
		double diff_L = blhs_new[i].L_Rad_ - blhs[i].L_Rad_;
		double diff_H = blhs_new[i].H_ - blhs[i].H_;
		double diff_B_deg, diff_L_deg;
		Rad2Deg(diff_B, diff_B_deg); Rad2Deg(diff_L, diff_L_deg);
		cout << setprecision(15) << diff_B_deg << "  " << diff_L_deg << "  " << diff_H << endl; 
	}
	cout << endl;
	for(auto blh : blhs_new) {
		cout << blh.B_Dms_.Deg_ << " " << blh.B_Dms_.Min_ << " " << blh.B_Dms_.Sec_ << "   " <<
				blh.L_Dms_.Deg_ << " " << blh.L_Dms_.Min_ << " " << blh.L_Dms_.Sec_ << "   " <<
				blh.H_ << endl;
	}
	cout << endl;
	// Gaussian Projection test
	GaussProj proj(config);
	proj.Proj_Batch();
	auto points = proj.GetPointsPlain();
	for (auto point : points)
		cout << fixed << setprecision(4) << point.X_ << "   " << point.Y_ << "   " << point.Z_ << endl;

	auto pts_blh = config.GetPoints();
	auto pts_blh_back = proj.BackProj_Batch();

	for (int i = 0; i < pts_blh.size(); ++i)
		pts_blh_back[i].H_ = pts_blh[i].H_;
	coor.SetPointsBLH(pts_blh_back);
	auto back = coor.BLH2XYZ_Batch();

	cout << endl << endl;
	for (int i = 0; i < back.size(); ++i) 
		cout << back[i].X_ - data[i].X_ << "   " <<
				back[i].Y_ - data[i].Y_ << "   " <<
				back[i].Z_ - data[i].Z_ << "   " << endl;
}