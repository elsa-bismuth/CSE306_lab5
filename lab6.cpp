#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <chrono>

#include "stb_image_write.h"

#include "fluid.h"
#include "edge.h"
#include "voronoi.h"


int main(){

    // voronoi1 : 24 points 
    // voronoi2 : 256 points
    // voronoi3 : 1024 points

    std::vector<Vector> points(1024); 

    for (int i=0; i < points.size(); i++){
        points[i][0] = rand()/ (double)RAND_MAX;
        points[i][1] = rand()/ (double)RAND_MAX;
        points[i][2] = 0;
    }

    VoronoiDiagram voronoi(points);
    voronoi.compute();
    voronoi.save("voronoi3.svg");

    return 0;
}