#ifndef polygon_h
#define polygon_h

#include <vector>
#include "vector.h"

class Polygon {  
public:
    std::vector<Vector> vertices;

    // Lab 7
    double area() const {
        if (vertices.size() < 3) {
            return 0; 
        }

        double result = 0;
        for (int i = 0; i < vertices.size(); i++) {
            const Vector& A = vertices[i];
            const Vector& B = vertices[(i + 1) % vertices.size()];
            result += A[0] * B[1] - A[1] * B[0];
        }

        return std::abs(result / 2);
    }

    // Lab 7
    double integrate_squared_distance(const Vector& P) const {
        if (vertices.size() < 3) {
            return 0;
        }

        double value = 0;
        for (int i = 1; i < vertices.size() - 1; i++) {
            Vector triangle[3] = {vertices[0], vertices[i], vertices[i + 1]};
            double local_value = 0;
            for (int k = 0; k < 3; k++) {
                for (int l = k; l < 3; l++) {
                    local_value += dot(triangle[k] - P, triangle[l] - P);
                }
            }

            Vector e1 = triangle[1] - triangle[0];
            Vector e2 = triangle[2] - triangle[0];
            double area_triangle = 0.5 * std::abs(e1[1] * e2[0] - e1[0] * e2[1]);
            value += local_value * area_triangle / 6.;
        }

        return value;
    }

    // Lab 8
    Vector centroid() const {
        Vector c(0, 0, 0);

        int N = vertices.size();
        double a = area();

        for (int i = 0; i < N; i++) {
            const Vector& currentVertex = vertices[i];
            const Vector& nextVertex = vertices[(i + 1) % N];

            c[0] += (currentVertex[0] + nextVertex[0]) * (currentVertex[0] * nextVertex[1] - nextVertex[0] * currentVertex[1]);
            c[1] += (currentVertex[1] + nextVertex[1]) * (currentVertex[0] * nextVertex[1] - nextVertex[0] * currentVertex[1]);
            c[2] += 0; // Z-coordinate not used in 2D polygons
        }

        return c / (6 * a);
    }

};  
 
#endif