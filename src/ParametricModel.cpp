#include "ParametricModel.h"
#include <algorithm>
#include <limits>

// Compute corrected a_hat and bx_hat using Eq. (23)
void ParametricModel::corrected_split_and_bx(const Params &p, double &a_hat, double &bx_hat) {
    double a = p[0], bx = p[1], by = p[2];
    double m = p[3], n = p[4];
    double d1 = p[5], d2 = p[6], d3 = p[7];

    // If all shifts equal the corrected split is trivial
    if (std::fabs(d1-d2)<1e-12 && std::fabs(d2-d3)<1e-12) {
        a_hat = a; bx_hat = bx; return;
    }

    double numA = a*std::cos(d2-d3)*std::cos(d1-d3) - bx*std::sin(d2-d3)*std::sin(d1-d3);
    double denom = std::cos(d1-d3)*std::cos(d2-d3) + std::sin(d1-d3)*std::sin(d2-d3);

    if (std::fabs(denom)<1e-12) { a_hat=a; bx_hat=bx; return; }

    a_hat = numA/denom;
    bx_hat = (a*std::sin(d1-d3)+bx*std::cos(d1-d3))/denom;
}

// Evaluate the parametric form
void ParametricModel::eval_alpha(double alpha, const Params &p, double &xout, double &yout) {
    double a = p[0], bx = p[1], by = p[2];
    double m = p[3], n = p[4];
    double d1 = p[5], d2 = p[6], d3 = p[7];
    double Hscale = p[8], Mscale = p[9];

    double a_hat, bx_hat;
    corrected_split_and_bx(p, a_hat, bx_hat);

    double x = a_hat * std::pow(std::cos(alpha+d1), m) + bx_hat * std::pow(std::sin(alpha+d2), n);
    double y = by * std::sin(alpha+d3);

    xout = Hscale*x;
    yout = Mscale*y;
}

// Sample curve uniformly at alpha across [0, 2PI]
void ParametricModel::sample_curve(const Params &p, size_t N, std::vector<double> &H, 
                                    std::vector<double> &M) {
    H.clear(); M.clear();
    H.reserve(N); M.reserve(N);
    const double TWO_PI=2.0* M_PI;
    for (size_t i=0; i < N ; i++) {
        double alpha = TWO_PI*(double(i)/double(N));
        double x,y;
        eval_alpha(alpha,p,x,y);
        H.push_back(x);
        M.push_back(y);
    }
}

// Loop area by numeric integration
double ParametricModel::loop_area(const Params &p, size_t N) {
    double area = 0.0;
    double TWO_PI = 2.0*M_PI;
    double delta_alpha = TWO_PI/double(N);

    double x0,y0;
    eval_alpha(0.0,p,x0,y0);

    for (size_t i=1 ; i <= N ; i++) {
        double alpha = delta_alpha*double(i);
        double x,y;
        eval_alpha(alpha,p,x,y);
        double dx = x-x0;
        double dy = y-y0;
        area += 0.5*(x0*dy - y0*dx);
        x0=x;
        y0=y;
    }
    return std::fabs(area);

}

// EvaluateM(H) at query H by sampling and interpolation
double ParametricModel::evaluate(double Hq, const Params &p, size_t N) {
    
    std::vector<double> H, M;
    sample_curve(p, N, H, M);

    // Nearest bracketing points
    for (size_t i = 1; i < H.size(); i++) {
        if((H[i-1] <= Hq && Hq <= H[i]) ||
            (H[i] <= Hq && Hq <= H[i-1])) {
                double t = (Hq - H[i-1]) / (H[i] - H[i-1]);
                return M[i-1] + t * (M[i] - M[i-1]);
            }
    }

    // Out of range: return the nearest endpoint
    if (Hq < H.front()) return M.front();
    if (Hq > H.back()) return M.back();
    return 0.0;
}