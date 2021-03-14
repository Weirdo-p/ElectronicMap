/*---------------------------------------------
   dijkscommon.h
   create on 13 Mar 2021 ZHUOXU WHU
---------------------------------------------*/
#ifndef _DIJKSCOMMON_H_
#define _DIJKSCOMMON_H_

 // std
#include <iostream>
#include <vector>
#include <limits>

// eigen
#include <Eigen/Core>

#define STATUS_NORMAL   0
#define STATUS_NO_ROUTE 1
#define STATUS_ERROR    2

using namespace std;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatNNd;
typedef vector<int>  VecInt;

#endif // _DIJKSCOMMON_H_