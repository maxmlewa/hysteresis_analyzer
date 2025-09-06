#ifndef LM_OPTIMIZER_H
#define LM_OPTIMIZER_H

#include <vector>
#include <functional>

// Options for the Levenberg-Marquardt optimization
struct LMOptions {
    int max_iters;
    double tol_grad;
    double tol_step;
    double lambda0;
    double lambda_up;
    double lambda_down;
    double fd_eps;

    LMOptions()
        : max_iters(300), tol_grad(1e-8), tol_step(1e-8), lambda0(1e-3),
          lambda_up(10.0), lambda_down(0.1), fd_eps(1e-6) {}
};


// Perform the LM Optimization on the given residual function
bool levenberg_marquardt(std::vector<double> &params,
                        const std::function<void(const std::vector<double>&,
                                            std::vector<double>&)> &residual_fn,
                        double &out_sse,
                        const LMOptions &opts = LMOptions());

#endif