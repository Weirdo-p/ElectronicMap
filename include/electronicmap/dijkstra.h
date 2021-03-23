/*---------------------------------------------
   dijkstra.h
   create on 09 Mar 2021 ZHUOXU WHU
---------------------------------------------*/

#ifndef _DIJKSTRA_H_
#define _DIJKSTRA_H_

#include "electronicmap/dijksconfig.h"
#include "electronicmap/dijkscommon.h"
#include <fstream>

class CDijkstra {   // dijkstra algorithm
public:
    CDijkstra();
    CDijkstra(CDijksConfig config);
    ~CDijkstra();

public: // main algoritm
    /********************************************************
     * find the nearest path from start vertex to end vertex
     * @param   verNumStart [in]    start vertex number
     * @param   verNumEnd   [in]    end vertex number
     * @return  route from start to end
    ********************************************************/
    VecInt FindPathNear(int verNumStart, int verNumEnd);

private: // helpers
    /***************************************************
     * search attach vertex
     * @param   vernum      [in]    vertex number
     * @param   startNum    [in]    start vertex number 
     * @return  vertex attached with vernum
    ***************************************************/
    VecInt SearchAttach(int vernum, int startNum);

    /********************************************************
     * update variables
     * @param   vertex  [in]    vertex index
     * @param   attach  [in]    index of attached with vertex
     * @param   mindist [in]    min distance for now 
     * @return
    ********************************************************/
    void UpdateState(int vertex, VecInt attach, double mindist);

    /***************************************************
     * search min distance in dist_ and isFind == false
     * @return index of min distance vertex
    ***************************************************/
    int DistMin();

    /*********************************************************
     * generate min dist from verNumStart to any other vertex
     * @param   verNumStart     [in]    number of start vertex
     * @return  STATUS_CODE
    *********************************************************/
    int GenerateMap(int verNumStart);

private:
    MatNNd  adj_matrix_;     /* Adjacency Matrix */
    int     vertex_num_;     /* number of vertex */
    bool*   isFind_;         /* set 1 if find the nearest
                                path to the origin vertex */
    double* dist_;           /* the nearest distance to
                                the origin vertex */
    int*    path_;           /* last vertex */
    double  mindist_;        /* min distance for now */
};

#endif // _DIJKSTRA_H_