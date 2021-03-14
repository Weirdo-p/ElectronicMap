#include <iostream>
using namespace std;

#include <Eigen/Core>
#include "electronicmap/dijksconfig.h"
#include "electronicmap/dijkstra.h"

int main(int argv, char** argc) {
    if(argv != 2) {
        cout << "usage:  ./dijkstest [filename]" << endl;
        return 0;
    }
    char* path = argc[1];
    remove("./data/result.txt"); remove("./data/path.txt");

    CDijksConfig config(path);
    CDijkstra runner(config);
    for(int i = 0; i < config.GetVerNum(); ++i)
        runner.FindPathNear(0, i);
    return 0;
}