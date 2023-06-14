#ifndef power_diagram_h
#define power_diagram_h

#include "polygon.h"

// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none") {
    FILE* f = fopen(filename.c_str(), "w+"); 
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    for (int i=0; i<polygons.size(); i++) {
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \""); 
        for (int j = 0; j < polygons[i].vertices.size(); j++) {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }
    fprintf(f, "</svg>\n");
    fclose(f);
}

class PowerDiagram {
public:

    Polygon disk;

    std::vector<Vector> points; 
    std::vector<double> weights;
    std::vector<Polygon> diagram; 

    PowerDiagram() {};
    PowerDiagram(std::vector<Vector>& pts, std::vector<double> &weights) {
        points = pts;
        weights = weights;

        // Lab 8
        const int N_disk = 50;
        disk.vertices.resize(N_disk);

        for (int i = 0; i < N_disk; i++){
            double theta = 2 * M_PI * i /(double)N_disk;
            disk.vertices[i][0]= cos(theta);
            disk.vertices[i][1]= sin(theta);
            disk.vertices[i][2]= 0;
        }
    }
    
    // Sutherland-Hodgman algorithm
    Polygon clip_polygon_by_bissector(const Polygon &poly, int index_0, int index_i, const Vector& P0, const Vector& Pi) {
        
        Polygon result;

        Vector M = (P0 + Pi) * 0.5; 
        Vector Mprime = M + ((weights[index_0]-weights[index_i])/(2.*(P0-Pi).norm2()))*(Pi-P0);

        for (int i=0; i < poly.vertices.size(); i++){
            const Vector &A = (i==0)? poly.vertices[poly.vertices.size() - 1]:poly.vertices[i-1];
            const Vector &B = poly.vertices[i];
            double t = dot(Mprime-A, Pi - P0) / dot(B - A, Pi - P0);
            Vector P = A + t*(B-A);

            // updated formula with p101
            if ((B-P0).norm2() - weights[index_0] < (B-Pi).norm2() - weights[index_i]){ // B is inside
                if ((A-P0).norm2() - weights[index_0] > (A-Pi).norm2() - weights[index_i] ){ // A is outside
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            }
            else {
                if ((A-P0).norm2() - weights[index_0]< (A-Pi).norm2() - weights[index_i]){ // A is inside
                    result.vertices.push_back(P);
                }
            }
        }
        return result;
    }

    // Lab 8
    Polygon clip_by_edge(const Polygon &poly, const Vector &u, const Vector &v) const {
       
        Polygon result;
        result.vertices.reserve(poly.vertices.size() + 1);
        Vector N(v[1] - u[1], -v[0] + u[0], 0.);

        for (int i = 0; i < poly.vertices.size(); i++) {
            const Vector &A = (i == 0) ? poly.vertices[poly.vertices.size() - 1] : poly.vertices[i - 1];
            const Vector &B = poly.vertices[i];
            double t = dot(u - A, N) / dot(B - A, N);
            Vector P = A + t * (B - A);

            if (dot(B - u, N) < 0) { // B is inside
                if (dot(A - u, N) > 0) { // A is outside
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            }
            else { // B is outside
                if ((dot(A - u, N) < 0)) { // A is inside
                    result.vertices.push_back(P);
                }
            }
        }
        return result;
    }

    // Lab 8
    Polygon intersect_with_disk(const Polygon& polygon, const Vector& center, double radius ) const {
        
        Polygon result(polygon);

        for (int i = 0; i < disk.vertices.size(); i++){
            const Vector &A = disk.vertices[i]*radius + center;
            const Vector &B = disk.vertices[(i+1)%disk.vertices.size()]*radius + center;
            result = clip_by_edge(result, A, B);
        }

        return result;
    }

    Polygon compute_diagram_cell(int idx){

        Polygon result;
        result.vertices.resize(4);
        result.vertices[0] = Vector(0, 0, 0);
        result.vertices[1] = Vector(0, 1, 0);
        result.vertices[2] = Vector(1, 1, 0);
        result.vertices[3] = Vector(1, 0, 0);

        for (int i=0; i<points.size(); i++){
            if (i==idx) {
                continue;
            }
            result = clip_polygon_by_bissector(result, idx, i, points[idx], points[i]);
        }

        // Lab 8
        // result = intersect_with_disk(result, points[idx], sqrt(weights[idx] - weights.back()));

        return result;

    }

    void compute() {

        diagram.resize(points.size());
        for (int i=0; i<points.size(); i++){
            diagram[i] = compute_diagram_cell(i);
        }
    }

    void save(std::string filename){
        save_svg(diagram, filename, "blue");
    }

};

#endif