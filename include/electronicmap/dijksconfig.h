/*---------------------------------------------
   dijksconfig.h
   create on 09 Mar 2021 ZHUOXU WHU
---------------------------------------------*/
#ifndef _DIJKSCONFIG_H_
#define _DIJKSCONFIG_H_

#include "electronicmap/dijkscommon.h"
#include <string>

class CDijksConfig {
public:
    CDijksConfig();
    CDijksConfig(char* filepath);
    CDijksConfig(MatNNd matrix);

public: // get function
    MatNNd GetAdjMat();
    int    GetVerNum();

private:
    void ReadFile(char* filepath);
    void InitAdjMatrix(string line);
    void SetAdjMatrix(string line);
    
private:
    MatNNd adj_matrix_;     /* Adjacency Matrix */
    int    vertex_num_;     /* number of vertex */
};
#endif // _DIJKSCONFIG_H_