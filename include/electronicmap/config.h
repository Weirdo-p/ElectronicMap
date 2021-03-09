/*---------------------------------------------
   config.h
   create on 09 Mar 2021 ZHUOXU WHU
---------------------------------------------*/
#ifndef _CONFIG_H_
#define _CONFIG_H_

#include "electronicmap/common.h"
#include <vector>

class CConfigCoors {
public:
    CConfigCoors();
    CConfigCoors(char* FilePath, DMS CentralLongi);

public: // set function
    void SetCentral(DMS central);
    void SetEllipsoid(EllipsoidType);

public: // get function
    DMS GetCentralDms();
    double GetCentralRad();
    vector<BLH> GetPoints();
    EllipsoidType GetEllipsoidType();

private:
    bool ReadFile(char* FilePath);

private:
    std::vector<BLH>    points_;        // undetermined points in BLH
    DMS                 central_dms_;   // central longitude in DMS
    double              central_rad_;   // central longitude in rad
    EllipsoidType       elli_;          // ellipsoid type
};

#endif // _CONFIG_H_