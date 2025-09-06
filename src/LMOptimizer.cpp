#include "LMOptimizer.h"
#include <cmath>
#include <limits>
#include <iostream>

/**
 * @brief Solve a linear system Ax = b using naive Gaussian elimination.
 * 
 * Sufficient for a small parameter count <= 20
 * 
 * @param A Coefficient matrix (nxn)
 * @param b RHS (n)
 * @param x Solution Vector (n)
 * @return true if solved successfully, false if singular
 */
//
static bool solve_linear(std::vector<std::vector<double>> A, std::vector<double> b,
                        std::vector<double> &x) {

    int n = static_cast<int>(b.size());
    x.assign(n, 0.0);

    // Forward elimination
    for (int i = 0 ; i < n ; i++) {
        // Pivoting
        int piv = i;
        double maxv = std::fabs(A[i][i]);

        for (int r = i+1; r < n; r++) {
            if (std::fabs(A[r][i]) > maxv) {
                maxv=std::fabs(A[r][i]);
                piv = r;
            }
        }
            
        if (piv != i) { 
            std::swap(A[i],A[piv]); 
            std::swap(b[i],b[piv]); 
        }

        double diag=A[i][i];

        if (std::fabs(diag)<1e-15) return false; // singular or not well formed

        // Eliminate below
        for (int r=i+1; r < n; r++) {
            double fac = A[r][i]/diag;

            for (int c = i; c < n; c++) {
                A[r][c] -= fac*A[i][c];
            }

            b[r] -= fac*b[i];
        }
    }

    // Back substitution
    for (int i = n-1; i >= 0; i--) {
        double s = b[i];

        for (int c = i + 1; c < n; c++) s -= A[i][c] * x[c];

        if (std::fabs(A[i][i]) < 1e-15) return false;

        x[i] = s/A[i][i];
    }

    return true;
}

/**
 * @brief Levenberg-Marguardt nonlinear least squares optimizer.
 * 
 * Attempts to minimize the sum of squared residuals by iteratively updating parameters.
 * Uses numerical Jacobian and a damped Gauss-Newton step.
 * 
 * @param params      [in/out] parameter vector to be optimized
 * @param residual_fn Function that computes residuals given params
 * @param out_see     [out] final sum of squared errors
 * @param opts        LM options (damping, max iterations, )
 * @return true if iterations complete, false if numerical failure
 */
bool levenberg_marquardt(std::vector<double> &params,
                         const std::function<void(const std::vector<double>&,
                                                  std::vector<double>&)> &residual_fn,
                         double &out_sse,
                         const LMOptions &opts) {

    size_t npar=params.size();

    // Evaluate initial residuals and SSE
    std::vector<double> res;
    residual_fn(params,res);
    size_t nres=res.size();
    double sse=0.0;

    for (double r:res) {
        sse+=r*r;
    }

    out_sse=sse;

    // Initialize the damping factor
    double lambda=opts.lambda0;

    for (int iter = 0; iter < opts.max_iters; iter++) {

        // Numerical Jacobian: J(i,j) = delta_r_i/delta_p_j approximated by finite differences
        std::vector<std::vector<double>> J(nres,std::vector<double>(npar, 0.0));

        for (size_t j = 0; j < npar; j++) {
            std::vector<double> p2 = params;
            double eps = opts.fd_eps * std::max(1.0, std::fabs(params[j]));
            p2[j] += eps;

            std::vector<double> r2;
            residual_fn(p2,r2);

            for (size_t i = 0; i < nres; i++) {

                J[i][j] = (r2[i] - res[i]) / eps;
            }
        }

        // Build normal equations: JTJ * dp = -JTr
        std::vector<std::vector<double>> JTJ(npar, std::vector<double>(npar, 0.0));
        std::vector<double> JTr(npar, 0.0);

        for (size_t i = 0; i < nres; i++) {

            for (size_t j = 0; j < npar; j++) {

                JTr[j] += J[i][j] * res[i];

                for (size_t k = 0; k < npar; k++) {
                    JTJ[j][k] += J[i][j] * J[i][k];
                }
            }
        }

        // Damping adjustment: JTJ + lambda*I (diagonal adjustment)
        for (size_t d = 0; d < npar; d++) JTJ[d][d] *= (1.0 + lambda);

        // RHS vector
        std::vector<double> rhs(npar);
        for (size_t j = 0; j < npar; j++) rhs[j] =- JTr[j];

        // Solve the linear system for update step
        std::vector<double> dp;
        if (!solve_linear(JTJ, rhs, dp)) return false; // numerical failure

        // Evaluating candidate parameters
        std::vector<double> cand = params;

        for (size_t j = 0; j < npar; j++) cand[j] += dp[j];

        std::vector<double> res_cand;
        residual_fn(cand,res_cand);

        double sse_cand=0.0;
        for (double r: res_cand) sse_cand += r*r;

        // Accept/Reject step
        if (sse_cand < sse) { // improvement, thus accepting update
            params = cand; 
            res = res_cand; 
            sse = sse_cand;
            lambda *= opts.lambda_down; // reducing the damping 

            if (lambda < 1e-20) lambda=1e-20;
        }

        else {
            // no improvement thus increasing damping
            lambda *= opts.lambda_up;
        }

        out_sse=sse;
    }
    
    return true;
}