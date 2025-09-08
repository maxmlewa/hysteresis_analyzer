#include "ModelFitter.h"
#include "LMOptimizer.h"
#include <cmath>
#include <limits>
#include <iostream>

// Implementation of ModelFitter

// Build residual vector for LM optimization.
void ModelFitter::build_residual(const std::vector<DataPoint>& data,
                                 const ParametricModel::Params &params,
                                 std::vector<double> &residuals) {
    residuals.resize(data.size());
    for (size_t i = 0; i < data.size(); i++) {
        residuals[i] = ParametricModel::evaluate(data[i].H, params) - data[i].M;
    }
}



/**
 * @brief Simple linear interpolation of M(H).
 *
 * Given two nearest sample points (H[i], M[i]) and (H[i+1], M[i+1]),
 * estimate M(Hq) by linear interpolation.
 */
double ModelFitter::interp_M_of_H(const std::vector<double>& Hs,
                                  const std::vector<double>& Ms,
                                  double Hq) {
    if (Hs.empty()) return 0.0;

    for (size_t i = 1; i < Hs.size(); i++) {
        if ((Hs[i-1] <= Hq && Hq <= Hs[i]) ||
            (Hs[i] <= Hq && Hq <= Hs[i-1])) {
            double t = (Hq - Hs[i-1]) / (Hs[i] - Hs[i-1]);
            return Ms[i-1] + t * (Ms[i] - Ms[i-1]);
        }
    }
    // Out of range → return nearest endpoint
    return (Hq < Hs.front()) ? Ms.front() : Ms.back();
}

/**
 * @brief Compute derived hysteresis metrics from a fitted model.
 *
 *   - Ms: value at very large |H|
 *   - Mr: magnetization at H=0
 *   - Hc: coercive field (H at which M=0)
 *   - Area: approximate integral over one loop branch
 */
void ModelFitter::compute_model_metrics(const ParametricModel::Params &p,
                                        FitResult &r) {
    // Saturation (evaluate at large field)
    r.Ms_model = ParametricModel::evaluate(1e6, p);

    // Remanence (at H=0)
    r.Mr_model = ParametricModel::evaluate(0.0, p);

    // Coercivity: brute-force scan for M=0 crossing
    double Hmin = -1e4, Hmax = 1e4;
    double step = (Hmax - Hmin) / 2000;
    double bestH = 0.0, bestDiff = std::numeric_limits<double>::max();

    for (double H = Hmin; H <= Hmax; H += step) {
        double M = ParametricModel::evaluate(H, p);
        if (std::fabs(M) < bestDiff) {
            bestDiff = std::fabs(M);
            bestH = H;
        }
    }
    r.Hc_model = bestH;

    // Loop area (energy loss proxy):
    // Integrate M(H) dH using trapezoidal rule.
    double area = 0.0;
    int N = 2000;
    double dH = 2 * Hmax / N;
    for (int i = 1; i < N; i++) {
        double H1 = -Hmax + (i-1) * dH;
        double H2 = -Hmax + i * dH;
        double M1 = ParametricModel::evaluate(H1, p);
        double M2 = ParametricModel::evaluate(H2, p);
        area += 0.5 * (M1 + M2) * (H2 - H1);
    }
    r.area_model = std::fabs(area);
}

/**
 * @brief Run Levenberg–Marquardt optimization with parameter mask.
 *
 * Mask lets us freeze certain parameters during staged fitting.
 *   mask[j] = true  → parameter j is optimized
 *   mask[j] = false → parameter j is fixed
 */
bool ModelFitter::run_lm_with_mask(const std::vector<DataPoint>& data,
                                   ParametricModel::Params &params,
                                   const std::vector<bool>& mask,
                                   double &sse) {
    // Build residual functor wrapper for LM
    auto residual_fn = [&](const std::vector<double>& p_vec,
                           std::vector<double>& residuals) {
        
        // Map masked parameter vector into full params
        ParametricModel::Params p = params;
        int idx = 0;
        for (size_t j = 0; j < mask.size(); j++) {
            if (mask[j]) {
                // update only free parameters
                p[j] = p_vec[idx++];
            }
        }
        build_residual(data, p, residuals);
    };

    // Flatten free parameters into vector
    std::vector<double> p_free;
    for (size_t j = 0; j < mask.size(); j++) {
        if (mask[j]) p_free.push_back(params[j]);
    }

    // Run LM optimizer
    LMOptions opts;
    opts.max_iters = 100;
    opts.lambda0   = 1e-3;
    opts.fd_eps    = 1e-6;

    bool ok = levenberg_marquardt(p_free, residual_fn, sse, opts);

    // Write back results
    if (ok) {
        int idx = 0;
        for (size_t j = 0; j < mask.size(); j++) {
            if (mask[j]) {
                params[j] = p_free[idx++];
            }
        }
    }
    return ok;
}

/**
 * @brief Main staged fitting routine.
 *
 * Strategy:
 *   Stage 1: Fit scale & slope params, hold phase shifts fixed
 *   Stage 2: Release phase shifts, refit
 *   Stage 3: Release all params, final global refinement
 *
 * This staged relaxation improves stability and avoids
 * local minima during nonlinear fitting.
 */
bool ModelFitter::fit_staged(const std::vector<DataPoint>& data,
                             ParametricModel::Params &initial_params,
                             FitResult &out) {
    ParametricModel::Params params = initial_params;
    double sse;

    // Stage 1: Optimize scale & slopes only
    {
        std::vector<bool> mask = {true, true, true, true, true,
                                  false, false, false, true, true};
        if (!run_lm_with_mask(data, params, mask, sse)) return false;
        std::cout << "[Stage 1] SSE=" << sse << std::endl;
    }

    // Stage 2: Release phase shifts (phi_x, phi_y, etc.)
    {
        std::vector<bool> mask = {true, true, true, true, true,
                                  true, true, true, true, true};
        if (!run_lm_with_mask(data, params, mask, sse)) return false;
        std::cout << "[Stage 2] SSE=" << sse << std::endl;
    }

    // Stage 3: Full refinement (all parameters free)
    {
        std::vector<bool> mask(10, true);
        if (!run_lm_with_mask(data, params, mask, sse)) return false;
        std::cout << "[Stage 3] SSE=" << sse << std::endl;
    }

    // Store final results
    out.final_params = params;
    out.sse = sse;
    compute_model_metrics(params, out);

    return true;
}