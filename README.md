# Project For Eletronic Map
This project is mainly to solve Dijkstra algorithm, and transfer a point(or a batch) from XYZ to BLH, or BLH to XYZ. And project them to Gaussian coordinate. 
## Data Format
### Coordinate Format
For a coordinate data file, we define symbol "#" as comments which at the beginning of a single line. The data format is as follow:

No. latitude-Deg latitude-Min latitude-Sec longitude-Deg longitude-Min longitude-Sec height

### Graph Format
It begins with a number which indicates how many vertex in this graph in a single line. And then, for each edges, we use "start-vertex end-vertex weight" to describe a graph.

## Compile
cd into your workfolder, and type these code:
```
    mkdir build
    cd build
    cmake..
    make
```
Then, you will see binary files in ${workfolder}/bin/, just run them!

