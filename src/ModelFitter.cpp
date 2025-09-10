// Robust, faster staged fitter for Lapshin parametric hysteresis model.
// - Downsamples the experimental data for fitting (keeps full data for metrics).
// - Performs staged LM (scale+shape -> exponents -> phase shifts -> refine all).
// - Uses conservative LM iteration limits to avoid long runs.
// - Emits SSE per stage for diagnostics.

#include "ModelFitter.h"
#include "LMOptimizer.h"
#include "ParametricModel.h"

#include <cmath>
#include <limits>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstddef>

// -----------------------------
// Utility: downsample data evenly
// -----------------------------
static std::vector<DataPoint> downsample_even(const std::vector<DataPoint>& data, size_t targetN) {
    if (data.size() <= targetN || targetN == 0) return data;
    std::vector<DataPoint> sampled;
    sampled.reserve(targetN);
    double step = static_cast<double>(data.size()) / static_cast<double>(targetN);
    double idx = 0.0;
    for (size_t k = 0; k < targetN; ++k) {
        size_t i = static_cast<size_t>(std::floor(idx));
        if (i >= data.size()) i = data.size() - 1;
        sampled.push_back(data[i]);
        idx += step;
    }
    return sampled;
}

// -----------------------------
// Basic residual builder
// r_i = M_model(H_i, params) - M_data_i
// Uses ParametricModel::evaluate(H, params) which interpolates model samples.
// -----------------------------
void ModelFitter::build_residual(const std::vector<DataPoint>& data,
                                 const ParametricModel::Params &params,
                                 std::vector<double> &residuals) {
    residuals.resize(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        double Hq = data[i].H;
        double Mmodel = ParametricModel::evaluate(Hq, params); // interpolation inside
        residuals[i] = Mmodel - data[i].M;
    }
}

// -----------------------------
// Linear interpolation helper (kept for potential use elsewhere)
// -----------------------------
double ModelFitter::interp_M_of_H(const std::vector<double>& Hs,
                                  const std::vector<double>& Ms,
                                  double Hq) {
    if (Hs.empty()) return 0.0;
    // find bracketing interval
    for (size_t i = 1; i < Hs.size(); ++i) {
        if ((Hs[i-1] <= Hq && Hq <= Hs[i]) || (Hs[i] <= Hq && Hq <= Hs[i-1])) {
            double denom = (Hs[i] - Hs[i-1]);
            if (std::fabs(denom) < 1e-15) return Ms[i];
            double t = (Hq - Hs[i-1]) / denom;
            return Ms[i-1] + t * (Ms[i] - Ms[i-1]);
        }
    }
    return (Hq < Hs.front()) ? Ms.front() : Ms.back();
}

// -----------------------------
// Compute derived metrics from fitted model
// -----------------------------
void ModelFitter::compute_model_metrics(const ParametricModel::Params &p,
                                        FitResult &r) {
    // Sample the parametric curve sufficiently fine
    std::vector<double> Hs, Ms;
    ParametricModel::sample_curve(p, 3000, Hs, Ms);

    // Ms_model = max |M|
    double maxabs = 0.0;
    for (double v : Ms) maxabs = std::max(maxabs, std::fabs(v));
    r.Ms_model = maxabs;

    // Mr_model (M at H ~ 0) via interpolation
    r.Mr_model = interp_M_of_H(Hs, Ms, 0.0);

    // Hc_model: find H where M crosses zero (closest to zero)
    double closest = std::numeric_limits<double>::max();
    double Hc = 0.0;
    for (size_t i = 1; i < Hs.size(); ++i) {
        if ((Ms[i-1] <= 0.0 && Ms[i] >= 0.0) || (Ms[i-1] >= 0.0 && Ms[i] <= 0.0)) {
            // linear root between i-1 and i
            double denom = (Ms[i] - Ms[i-1]);
            double t = (std::fabs(denom) < 1e-15) ? 0.0 : (0.0 - Ms[i-1]) / denom;
            double Hroot = Hs[i-1] + t * (Hs[i] - Hs[i-1]);
            if (std::fabs(Hroot) < closest) { closest = std::fabs(Hroot); Hc = Hroot; }
        }
    }
    r.Hc_model = Hc;

    // Area: use parametric loop area routine
    r.area_model = ParametricModel::loop_area(p, 3000);
}

// -----------------------------
// Helper: check params for numerical issues
// -----------------------------
static bool params_valid(const ParametricModel::Params &p) {
    for (double v : p) {
        if (!std::isfinite(v)) return false;
        // optional: clamp extremely large values
        if (std::fabs(v) > 1e12) return false;
    }
    return true;
}

// -----------------------------
// Run LM with a mask of free parameters.
// We limit LM iterations for speed and use the provided residual function.
// -----------------------------
bool ModelFitter::run_lm_with_mask(const std::vector<DataPoint>& data,
                                   ParametricModel::Params &params,
                                   const std::vector<bool>& mask,
                                   double &sse) {
    // Build vector of free parameters
    std::vector<double> p_free;
    p_free.reserve(params.size());
    for (size_t j = 0; j < params.size(); ++j) {
        if (mask[j]) p_free.push_back(params[j]);
    }

    // If no free parameters, just compute residuals and sse
    if (p_free.empty()) {
        std::vector<double> res;
        build_residual(data, params, res);
        double s = 0.0; for (double r : res) s += r*r;
        sse = s;
        return true;
    }

    // Prepare LM options: conservative defaults to avoid long runs
    LMOptions opts;
    opts.max_iters = 60;   // keep small so fit is quick; increase if needed
    opts.fd_eps = 1e-6;
    opts.lambda0 = 1e-3;
    opts.lambda_up = 10.0;
    opts.lambda_down = 0.1;

    // Residual function that maps p_free -> full params -> residuals
    auto residual_fn = [&](const std::vector<double>& pf, std::vector<double>& residuals) {
        ParametricModel::Params p_full = params;
        size_t idx = 0;
        for (size_t j = 0; j < mask.size(); ++j) {
            if (mask[j]) {
                if (idx < pf.size()) p_full[j] = pf[idx++];
            }
        }
        // Build residuals on provided data
        build_residual(data, p_full, residuals);
    };

    // Run LM
    bool ok = levenberg_marquardt(p_free, residual_fn, sse, opts);

    // Copy free back to params (even if ok==false, p_free may have been modified)
    size_t idx = 0;
    for (size_t j = 0; j < mask.size(); ++j) {
        if (mask[j]) {
            if (idx < p_free.size()) params[j] = p_free[idx++];
        }
    }

    // Validate numeric sanity
    if (!params_valid(params)) {
        std::cerr << "[ModelFitter] numeric instability in parameters detected\n";
        return false;
    }

    return ok;
}

// -----------------------------
// Main staged fitting routine
// -----------------------------
bool ModelFitter::fit_staged(const std::vector<DataPoint>& data,
                             ParametricModel::Params &initial_params,
                             FitResult &out) {
    if (data.empty() || initial_params.size() < 10) {
        std::cerr << "[ModelFitter] invalid input to fit_staged\n";
        return false;
    }

    // Downsample experimental data for the optimizer only (keep full data for final metrics)
    const size_t FIT_POINTS = 300; // target number of points used by LM
    std::vector<DataPoint> fit_data = downsample_even(data, FIT_POINTS);

    // Start from provided initial guess but make a better automatic scaling for Hscale & Mscale
    ParametricModel::Params params = initial_params;

    // Estimate Hscale and Mscale from data if they appear to be zero or tiny
    double maxH = 0.0, maxM = 0.0;
    for (const auto &d : data) { maxH = std::max(maxH, std::fabs(d.H)); maxM = std::max(maxM, std::fabs(d.M)); }
    if (maxH <= 0.0) maxH = 1.0;
    if (maxM <= 0.0) maxM = 1.0;

    // If user-provided initial_params don't set Hscale/Mscale meaningfully, set them
    if (params[8] == 0.0) params[8] = maxH;
    if (params[9] == 0.0) params[9] = maxM;

    // If a (param 0) is zero, give a modest starting value
    if (std::fabs(params[0]) < 1e-8) params[0] = 0.9;
    if (std::fabs(params[1]) < 1e-8) params[1] = 1.0;
    if (std::fabs(params[2]) < 1e-8) params[2] = 1.0;
    if (std::fabs(params[3]) < 1e-8) params[3] = 3.0;
    if (std::fabs(params[4]) < 1e-8) params[4] = 1.0;

    double sse = 0.0;
    bool ok;

    // ---------- Stage 1: fit Hscale, Mscale, and basic shape (a, bx, by) ----------
    {
        std::vector<bool> mask(params.size(), false);
        mask[0] = true; // a
        mask[1] = true; // bx
        mask[2] = true; // by
        mask[8] = true; // Hscale
        mask[9] = true; // Mscale

        std::cout << "[ModelFitter] Stage 1: fit scale + basic shape on " << fit_data.size() << " pts\n";
        ok = run_lm_with_mask(fit_data, params, mask, sse);
        std::cout << "  Stage 1 SSE = " << sse << "  (ok=" << (ok? "true":"false") << ")\n";
        if (!ok) {
            // continue attempt: proceed but warn
            std::cerr << "[ModelFitter] Stage 1 failed to converge cleanly; continuing to Stage 2\n";
        }
    }

    // ---------- Stage 2: fit exponents (m,n) ----------
    {
        std::vector<bool> mask(params.size(), false);
        mask[3] = true; // m
        mask[4] = true; // n
        // keep previously fit params free as well for joint re-fit stability
        mask[0] = true; mask[1] = true; mask[2] = true; mask[8] = true; mask[9] = true;

        std::cout << "[ModelFitter] Stage 2: fit exponents (m,n)\n";
        ok = run_lm_with_mask(fit_data, params, mask, sse);
        std::cout << "  Stage 2 SSE = " << sse << "  (ok=" << (ok? "true":"false") << ")\n";
        if (!ok) {
            std::cerr << "[ModelFitter] Stage 2 had issues; proceeding to Stage 3\n";
        }
    }

    // ---------- Stage 3: fit phase shifts (Δα1..Δα3) ----------
    {
        std::vector<bool> mask(params.size(), false);
        mask[5] = true; mask[6] = true; mask[7] = true; // Δα1..Δα3
        // keep shape + exponents + scales free to allow small adjustments
        mask[0] = mask[1] = mask[2] = mask[3] = mask[4] = true;
        mask[8] = mask[9] = true;

        std::cout << "[ModelFitter] Stage 3: fit phase shifts\n";
        ok = run_lm_with_mask(fit_data, params, mask, sse);
        std::cout << "  Stage 3 SSE = " << sse << "  (ok=" << (ok? "true":"false") << ")\n";
        if (!ok) {
            std::cerr << "[ModelFitter] Stage 3 issues; continuing to final refine\n";
        }
    }

    // ---------- Final: refine all parameters ----------
    {
        std::vector<bool> mask(params.size(), true);
        std::cout << "[ModelFitter] Final: refine all parameters (short run)\n";
        ok = run_lm_with_mask(fit_data, params, mask, sse);
        std::cout << "  Final SSE = " << sse << "  (ok=" << (ok? "true":"false") << ")\n";
        if (!ok) {
            std::cerr << "[ModelFitter] Final refinement failed\n";
            // We may still accept the current params if numerically valid
            if (!params_valid(params)) return false;
        }
    }

    // Save results and compute metrics on full parametric curve
    out.final_params = params;
    out.sse = sse;
    compute_model_metrics(params, out);

    return true;
}