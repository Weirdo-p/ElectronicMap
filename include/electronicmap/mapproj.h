/*---------------------------------------------
   mapproj.h
   create on 09 Mar 2021 ZHUOXU WHU
---------------------------------------------*/
#ifndef _MAPPROJ_H_
#define _MAPPROJ_H_

#include "electronicmap/common.h"
#include "electronicmap/config.h"
#include "electronicmap/coordinates.h"
#include <vector>
using namespace std;

class GaussProj {
public:
    GaussProj();
    GaussProj(CConfigCoors config);

public:
    /**********************************
     * to execute Gauss projection
     * result stores in points_plain_
     * @return  true if success 
    **********************************/
    XYZ Proj(BLH blh);

    /********************************************
     * to execute Gauss projection BATCH VERSION
     * result stores in points_plain_
     * @return  true if success 
    ********************************************/
    bool Proj_Batch();

    /**********************************************
     * Gaussian Back calculation
     * @param   xyz     Gaussian Plain coordinate
     * @return  BLH
    **********************************************/
    BLH BackProj(XYZ xyz);

    /**********************************************
     * Gaussian Back calculation Batch_Version
     * @return  BLH
    **********************************************/
    vector<BLH> BackProj_Batch();

    /******************************************
     * get central longitude
     * @param   longi     longitude
     * @return  central longitude
    ******************************************/
    int GetL0(DMS central);

private: // helpers for projection
    bool GetN(double B,double &N);

    /******************************************************
     * eta = ee2 * (cosB) ^ 2
     * @param   B   latitude of undetermined point(radians)
     * @return  eta
    ******************************************************/
    double Geteta(double B);

    /******************************
     * get coefficient
     * @return  5 * 1 double array
    ******************************/
    double* GetCoef();

    /************************************
     * @return  x coordinate in Gaussian
    ************************************/
    double Getx(double B, double t, double eta, double l,double N, double X);

    /************************************
     * @return  y coordinate in Gaussian
    ************************************/
    double Gety(double B, double t, double eta, double l,double N);

public: // helpers for back calculation
    /****************************************
     * @param   x   x coor of Gaussian plain
     * @return  底点维度
    ****************************************/
    double GetBf(double x);

    /***********************************************
     * to calculate radius of curvature of meridian
     * @param   B   底点维度
     * @return  radius of curvature of meridian
    ***********************************************/
    double GetM(double B);

    double GetL(double y, double Bf, double Nf, double tf, double etaf);
    double GetB(double Bf, double y, double tf, double Mf, double Nf, double etaf);
public: // get functions
    vector<XYZ> GetPointsPlain();


private:
    vector<BLH> points_blh_;
    vector<XYZ> points_xyz_;
    vector<XYZ> points_plain_;
    DMS         central_dms_;
    Ellipsoid   elli_;
    double      central_rad_;
    CCoors      coor;
};

#endif // _MAPPROJ_H_