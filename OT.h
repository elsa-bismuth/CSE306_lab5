#ifndef ot_h
#define ot_h

#include "power_diagram.h"
#include "lbfgs.c"

#define VOLUME_FLUID 3.0
#define VOLUME_AIR 1.0

class OT{
public:
    OT(const std::vector<Vector> &pts, const std::vector<double> &lambdas){
        this->pts = pts;
        this->lambdas=lambdas;
    };

    // from sample.cpp
    static lbfgsfloatval_t _evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        return reinterpret_cast<OT*>(instance)->evaluate(x, g, n, step);
    }

    lbfgsfloatval_t evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        lbfgsfloatval_t fx = 0.0;
     
        for (int i = 0; i < n ; i++) {
            solution.weights[i] = x[i];
        }
        solution.compute();

        double s1 = 0;
        double s2 = 0;
        double s3 = 0;
        double estimated_volume_fluid = 0;

        // p102-103
        for (int i = 0; i < n ; i++) {
            double cell_area = solution.diagram[i].area();
            //g[i]= -(this->lambdas[i] - cell_area);
            g[i] = static_cast<lbfgsfloatval_t>(-(this->lambdas[i] - cell_area));
            s3 += this->lambdas[i]*x[i];
            s2 -= (x[i]*cell_area);
            s1 += solution.diagram[i].integrate_squared_distance(solution.points[i]);
            estimated_volume_fluid += cell_area;
        }
        
        fx  = s1 + s2 + s3;

        // Lab 8
        double estimated_volume_air = 1. - estimated_volume_fluid;
        fx += x[n-1] * (VOLUME_AIR - estimated_volume_air);
        g[n-1] = -(VOLUME_AIR - estimated_volume_air);

        return -fx;

    }

    static int _progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        return reinterpret_cast<OT*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
    }

    int progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        for (int i = 0;i < n;i++) {
            solution.weights[i] = x[i];
        }
        solution.compute();

        double max_diff=0;
        for(int i=0; i< n; i++){
            double current_area = solution.diagram[i].area();
            double desired_area = lambdas[i];
            max_diff = std::max(max_diff, std::abs(current_area-desired_area));
        }

        std::cout<<"fx:"<< fx<<"\tmax difference = "<< max_diff<<"\t gnorm"<< gnorm <<std::endl;

        return 0;
    }

    // solves the optimal transport problem
    void solve(){
        solution.points = this->pts;
        solution.weights.resize(this->lambdas.size()); 
        std::fill(solution.weights.begin(), solution.weights.end(), this->lambdas);

        double fx = 0;

        // LBFGS function from sample.cpp 
        int ret = lbfgs(this->pts.size(), &solution.weights[0], &fx, _evaluate, _progress, this, NULL);
        solution.compute();

    }

    std::vector<Vector> pts;
    std::vector<double> lambdas;
    PowerDiagram solution;

};

#endif