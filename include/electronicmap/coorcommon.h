/*---------------------------------------------
   coorcommon.h
   create on 09 Mar 2021 ZHUOXU WHU
---------------------------------------------*/
#ifndef _COORCOMMON_H_
#define _COORCOMMON_H_

#include <iostream>
#include <math.h>

#define PI (atan(1) * 4)

enum EllipsoidType {CGCS2000, WGS84};   /* support WGS84 CGCS2000 */

struct DMS {                /* degree form */
public:
    int     Deg_;           /* degree */
    int     Min_;           /* minute */
    double  Sec_;           /* second */

public:
    DMS();
    DMS(int deg, int min, double sec);
private:
    void DMSCheck();
};

struct Giant {              /* big number type */
public:
    int     Integer_;       /* integer part */
    double  Frac_;          /* fractional part */

public:
    Giant();
    Giant(int inte, double frac);
};

struct BLH {                /* blh type */
public:
    double  B_Rad_;         /* latitude in rad */
    double  L_Rad_;         /* longitude in rad */
    double  H_;             /* geodetic height in meter*/
    DMS     B_Dms_;         /* latitude in DMS */
    DMS     L_Dms_;         /* longitude in DMS */

public:
    BLH();
    BLH(double B, double L, double H);
    BLH(DMS B, DMS L, double H);
};

struct XYZ {                /* ecef coordinates */
public:
    double  X_;              /* X coordinate */
    double  Y_;              /* Y coordinate */
    double  Z_;              /* Z coordinate */

public:
    XYZ();
    XYZ(double x, double y, double z);
};

struct Ellipsoid {          /* ellipsoid type(default WGS84) */
public:
    double a_;              /* major semi axis */
    double b_;              /* minor semi axis */
    double c_;              /* helper param */
    double e2_;             /* square of first eccentricity */
    double alpha_ ;         /* oblateness */
    double eprime2_;        /* second eccentricity */
    EllipsoidType type_;    /* ellipsoid type */

public:
    Ellipsoid();
    Ellipsoid(EllipsoidType type);
    Ellipsoid(double a, double b);

public:
    /**********************************
     * set ellipsoid params
     * @param   type    ellipsoid type
     * @return  true if success
    **********************************/
    bool SetEllipsoidParam(EllipsoidType type);

private:
    /**********************************
     * calculate ellipsoid params
     * @param   a   major semi axis
     * @param   b   minor semi axis
     * @return true if success
    **********************************/
    bool CalculateParam(double a, double b);
};

/***********************************
 * transform degree to radians
 * @param   deg     [In]    degree
 * @param   rad     [out]   radians
 * @return 
***********************************/
void Deg2Rad(double deg, double &rad);

/*****************************************
 * transform degree type from dms to deg
 * @param   dms     [in]    degree in dms
 * @param   deg     [out]   degree in deg
*****************************************/
void DMS2Deg(DMS dms, double &deg);

/***********************************
 * transform radians to degree
 * @param   rad     [In]    degree
 * @param   deg     [out]   radians
 * @return 
***********************************/
void Rad2Deg(double rad, double& deg);

/***********************************
 * transform rad to dms
 * @param   rad     [in]    radians
 * @param   dms     [out]   DMS
 * @return 
***********************************/
void Rad2DMS(double rad, DMS& dms);

/********************************
 * latitude validity check
 * @param   B   [in]    latitude
 * @return  true if legal
********************************/
bool CheckB(DMS B);

/********************************
 * longitude validity check
 * @param   L   [in]    longitude
 * @return  true if legal
********************************/
bool CheckL(DMS L);


#endif // _COMMON_H_