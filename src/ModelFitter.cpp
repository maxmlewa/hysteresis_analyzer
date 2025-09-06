#include "ModelFitter.h"
#include <cmath>
#include <limits>

// Modell Function
double ModelFitter::modelFunction(double H, double Ms, double a, double b, double c, double d) {
    
    return Ms * ( a * tanh(b * H) + (c * H)/ (1.0 + d * H * H));

}

// Simple gradient descent optimizer for the nonlinear fit
// to be optimized further

std::array<double, 4> ModelFitter::fitModel(const std::vector<DataPoint>& data, double Ms) {
    
    std::array<double, 4> params = {1.0, 1.0, 0.0, 0.0}; // initial guess
    double lr = 1e-3; // learning rate
    size_t maxIter = 500;

    for (size_t iter = 0; iter < maxIter; iter++) {
        std::array<double, 4> grad = {0, 0, 0, 0};
        double error = 0.0;

        for (size_t i = 0; i< data.size(); i++){
            double pred = modelFunction(data[i].H, Ms, params[0], params[1], params[2], params[3]);
            double diff = pred - data[i].M;
            error += diff *diff;

            // finite difference gradient
            double eps = 1e-6;
            for (int j = 0; j < 4; j++) {
                std::array<double, 4> temp = params;
                temp[j] += eps;
                double pred2 = modelFunction(data[i].H, Ms, temp[0], temp[1], temp[2], temp[3]);
                double grad_j = (pred2 - pred)/eps * diff;
                grad[j] += grad_j;
            }
        }

        for (int j = 0; j < 4; j++) {
            params[j] -= lr*grad[j];
        }
    }
}