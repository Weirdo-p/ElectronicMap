#ifndef _COORS_H_
#define _COORS_H_

#include "electronicmap/common.h"
#include "electronicmap/config.h"
#include <vector>
using namespace std;

class CCoors {
public:
    CCoors();
    CCoors(CConfigCoors config);

public:
    void BLH2XYZ();             // transform blh to xyz

private:
    vector<BLH> points_blh_;    // points in geodetic coordinate
    vector<XYZ> points_xyz_;    // points in ecef coordinate
    Ellipsoid   elli_;          // ellipsoid params
};

#endif // _COORS_H_