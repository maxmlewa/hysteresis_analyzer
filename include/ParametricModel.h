#ifndef PARAMETRIC_MODEL_H
#define PARAMETRIC_MODEL_H

#include <vector>
#include <array>
#include <cmath>
#ifndef M_PI;
#define M_PI 3.14159265358979323846
#endif

// Lapshin's improved parametric model (Rev. Sci. Instrum. 91, 065106 (2020))
class ParametricModel {
public:
    // Parameter vector layout:
    // [0] = a, [1] = bx, [2] = by, [3] = m, [4] = n,
    // [5] = delta_alpha1, [6] = delta_alpha2, [7] = delta_alpha3
    // [8] = Hscale, [9] = Mscale

    typedef std::vector<double> Params;

    // Evaluate the parametric equations x(alpha), y(alpha)
    static void eval_alpha(double alpha, const Params &p, double &xout, double &yout);

    // Sample the curve uniformly in alpha
    static void sample_curve(const Params &p, size_t N,
                            std::vector<double> &H, std::vector<double> &M);
            
    // Correct a_hat and b_hat per Equation (23) in paper
    static void corrected_split_and_bx(const Params &p, double &a_hat, double &bx_hat);

    // Loop area (numeric integration)
    static double loop_area(const Params &p, size_t N = 2000);
};

#endif