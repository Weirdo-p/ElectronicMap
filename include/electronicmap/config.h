#ifndef _CONFIG_H_
#define _CONFIG_H_

#include "electronicmap/common.h"
#include <vector>

class CConfigCoors {
public:
    CConfigCoors();
    CConfigCoors(char* FilePath, DMS CentralLongi);

private:
    std::vector<BLH>    points_;        // undetermined points in BLH
    DMS                 central_dms_;   // central longitude in DMS
    double              central_rad_;   // central longitude in rad
};

#endif // _CONFIG_H_