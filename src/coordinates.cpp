#include "electronicmap/coordinates.h"

CCoors::CCoors() { }

CCoors::CCoors(CConfigCoors config) {
    points_blh_ = config.GetPoints();
    elli_.SetEllipsoidParam(config.GetEllipsoidType());
}

