#ifndef POLYNOMIAL_MODEL_H
#define POLYNOMIAL_MODEL_H

#include "DataLoader.h"
#include <vector>
#include <string>

struct PolynomialFitResult {
    int degree;
    std::vector<double> coeffs_up;
    std::vector<double> coeffs_down;
    double Hmin, Hmax;
};

class PolynomialModel {
public:
    static PolynomialFitResult fit(const std::vector<DataPoint>& data, int degree = 9);
    static double eval(const std::vector<double>& coeffs, double H);

    static double computeSaturation(const PolynomialFitResult& r);
    static double computeRemanence(const PolynomialFitResult& r);
    static double computeCoercivity(const PolynomialFitResult& r);
    static double computeLoopArea(const PolynomialFitResult& r, int N = 2000);

    static bool write_csv(const std::string& path, const PolynomialFitResult& r, int N = 1000);
};

#endif