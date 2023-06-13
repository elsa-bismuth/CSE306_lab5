#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <chrono>
#include <map>
#include <sstream>
#include "stb_image_write.h"

#include "OT.h"
#include "edge.h"

int main() {

    // powerdiagram1 : lambdas[i] = 1
    // powerdiagram2 : rand()/(double)RAND_MAX
    // powerdiagram3 : 1./points.size()

    std::vector<Vector> points(256);
    std::vector<double> lambdas(256);

    Vector C = Vector(0.5,0.5,0.);
    double total = 0;

    for (int i = 0; i < points.size(); i++){
        points[i][0] = rand() / (double)RAND_MAX;
        points[i][1] = rand() / (double)RAND_MAX;
        points[i][2] = 0;
        //lambdas[i] = 1;
        //lambdas[i] = rand()/(double)RAND_MAX;
        Vector diff = C - points[i];
        lambdas[i] = std::exp(-pow(diff.norm(), 2.) / 0.02);
        total += lambdas[i];
    }

    // only for powerdiagram3 with semi-discrete optimal transport with LBFGS
    for (int i = 0; i < points.size(); i++){
        lambdas[i] /= total;
    }


    OT ot(points, lambdas);

    // Chrono
    auto start = std::chrono::high_resolution_clock::now();
    ot.solve();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout<<"duration = "<< duration.count() <<std::endl;

    ot.solution.save("powerdiagram3.svg");
    return 0;

}
