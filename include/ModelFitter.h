#ifndef MODEL_FITTER_H
#define MODEL_FITTER_H

#include "DataLoader.h"
#include "ParametricModel.h"
#include <vector>
#include <array>


// Container for the results of hysteresis fitting
struct FitResult {
    ParametricModel::Params final_params;
    double sse;

    // derived characteristics
    double Ms_model;
    double Mr_model;
    double Hc_model;
    double area_model;
};

// Provides staged nonlinear fitting of Lapshin's parametric hysteresis model
class ModelFitter {
public:

    /**
     * @brief Perform stged fitting of Lapshin's parametric model
     * 
     * @param data              Experimental hysteresis loop data
     * @param initial_params    Initial guess for the parameters
     * @param out               Struct filled with the parameters and metrics
     * @return true if fit converged, false on failure
     */
    static bool fit_staged(const std::vector<DataPoint>& data,
                           ParametricModel::Params &initial_params,
                           FitResult &out);
    
private:

        /**
         * @brief Construct residual vector for LM optimization
         * 
         * Residuals: r_i = M_model(H_i, params) - M_exp(H_i)
         * 
         * @param data          Experimantal data points
         * @param params        Model parameter set
         * @param residuals     Output vector of residuals
         */
        static void build_residual(const std::vector<DataPoint>& data,
                                   const ParametricModel::Params &params,
                                   std::vector<double> &residuals);

        /**
        * @brief Linear interpolation helper: M(H) from discrete samples.
        *
        * Useful for locating coercivity Hc (zero crossing of M).
        *
        * @param Hs   Sampled field values
        * @param Ms   Corresponding magnetization values
        * @param H    Query field
        * @return Interpolated magnetization
        */
        static double interp_M_of_H(const std::vector<double>& Hs,
                                    const std::vector<double>& Ms, double H);

        /**
        * @brief Compute derived loop metrics (Ms, Mr, Hc, area) from model params.
        *
        * This is run after fitting to summarize physical results
        * and enable comparison against simple analyzer outputs.
        *
        * @param p  Model parameters
        * @param r  Fit result struct (to be populated)
        */
        static void compute_model_metrics(const ParametricModel::Params &p,
                                          FitResult &r);


        /**
        * @brief Run LM optimization with a parameter mask.
        *
        * Mask allows staged fitting by fixing some parameters:
        *   - mask[j] = true  → optimize parameter j
        *   - mask[j] = false → hold parameter j fixed
        *
        * @param data    Experimental data
        * @param params  [in/out] parameters (updated in place)
        * @param mask    Which parameters are free to move
        * @param sse     [out] final sum of squared errors
        * @return true if LM converged, false otherwise
        */
        static bool run_lm_with_mask(const std::vector<DataPoint>& data,
                                     ParametricModel::Params &params,
                                     const std::vector<bool>& mask,
                                     double &sse);
};


#endif