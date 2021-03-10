/*---------------------------------------------
   coordinates.h
   create on 09 Mar 2021 ZHUOXU WHU
---------------------------------------------*/

#ifndef _COORS_H_
#define _COORS_H_

#include "electronicmap/common.h"
#include "electronicmap/config.h"
#include <vector>
using namespace std;

class CCoors { // coordinate related
public:
    CCoors();
    CCoors(CConfigCoors config);

public: // coordinate transformation
    /******************************************************
     * transform from blh to ecef coordinate(batch version)
     * @return  points in ecef coordinate,
     *          also store them in this->points_xyz_
    ******************************************************/
    vector<XYZ> BLH2XYZ_Batch();

    /******************************************************
     * transform from blh to ecef coordinate(batch version)
     * @return  points in ecef coordinate,
     *          also store them in this->points_xyz_
    ******************************************************/
    vector<BLH> XYZ2BLH_Batch();

    /******************************************************
     * transform a single point from blh to ecef coordinate
     * @param   blh [in]    geodetic coordinate
     * @param   xyz [out]   ecef coordinate
     * @return  point in ecef coordinate
    ******************************************************/
    bool BLH2XYZ(BLH blh, XYZ& xyz);

    /******************************************************
     * transform a single point from blh to ecef coordinate
     * @param   xyz [in]    ecef coordinate
     * @param   blh [in]    geodetic coordinate
     * @return  point in geodetic coordinate
    ******************************************************/
    bool XYZ2BLH(XYZ xyz, BLH &blh);

public: // get function
    vector<BLH> GetPointsBLH();
    vector<XYZ> GetPointsXYZ();

public: // set function
    void SetEllipsoid(EllipsoidType type);
    void SetPointsBLH(vector<BLH> points_blh);

private: // helper
    /****************************************************
     * 计算卯酉圈曲率半径
     * @param   B           [in]   latitude in radians
     * @param   ellipsoid   [in]   ellipsoid params
     * @param   N           [out]  卯酉圈曲率半径
     * @return  true if success
    ****************************************************/
    bool GetN(const double B, double &N);

private:
    vector<BLH> points_blh_;    // points in geodetic coordinate
    vector<XYZ> points_xyz_;    // points in ecef coordinate
    Ellipsoid   elli_;          // ellipsoid params
};

#endif // _COORS_H_