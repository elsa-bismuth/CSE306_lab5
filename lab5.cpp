#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <cstring>

std::random_device rd;
std::mt19937 engine(rd());
std::normal_distribution<dle> normal(0.0, 1.0);

void sliced_matching(const double* A, const double* B, int W, int H, double* result) {
    
    memcpy(result, A, W * H * 3 * sizeof(double));

    for (int iter = 0; iter < 100; iter++) {
        double xDir = normal(engine);
        double yDir = normal(engine);
        double zDir = normal(engine);
        double norm = sqrt(xDir * xDir + yDir * yDir + zDir * zDir);
        xDir /= norm;
        yDir /= norm;
        zDir /= norm;

        std::vector<std::pair<double, int> > resultSorted(W * H);

        std::vector<double> projA(W * H);
        std::vector<double> projB(W * H);
        for (int i = 0; i < W * H; i++) {
            projA[i] = A[i * 3] * xDir + A[i * 3 + 1] * yDir + A[i * 3 + 2] * zDir;
            projB[i] = B[i * 3] * xDir + B[i * 3 + 1] * yDir + B[i * 3 + 2] * zDir;
            resultSorted[i] = std::make_pair(projA[i], i);
        }

        std::sort(resultSorted.begin(), resultSorted.end());

        for (int i = 0; i < W * H; i++) {
            int index = resultSorted[i].second;
            double diff = projB[resultSorted[i].second] - resultSorted[i].first;
            result[index * 3] += diff * xDir;
            result[index * 3 + 1] += diff * yDir;
            result[index * 3 + 2] += diff * zDir;
        }
    }
}
