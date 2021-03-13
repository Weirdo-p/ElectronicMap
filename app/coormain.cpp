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
    char* path = "./data/points.txt";
    DMS central(113, 0, 0);
    CConfigCoors config(path, central);
    config.SetCentral(central);
    CCoors coor(config);
    auto data = coor.BLH2XYZ_Batch();

	// Gaussian Projection test
	GaussProj proj(config);
	proj.Proj_Batch();
	auto points = proj.GetPointsPlain();
	for (auto point : points)
		cout << setprecision(15) << point.X_ << "   " << point.Y_ << "   " << point.Z_ << endl;

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