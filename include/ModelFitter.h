#ifndef MODEL_FITTER_H
#define MODEL_FITTER_H

#include "DataLoader.h"
#include <vector>
#include <array>

class ModelFitter {
public:
    // Parameters: {a, b, c, d} for the model
    static std::array<double, 4> fitModel(const std::vector<DataPoint>& data,
                                            double Ms = 1.0);
    
    static double modelFunction(double H, double Ms, double a, double b, double c, double d);

};


#endif