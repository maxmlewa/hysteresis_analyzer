#include "HysteresisAnalyzer.h"
#include <cmath>
#include <limits>

double HysteresisAnalyzer::computeSaturation(const std::vector<DataPoint>& data) {
    double maxM = 0.0;
    for(size_t i = 0; i< data.size(); i++) {
        if(fabs(data[i].M > maxM)) maxM = fabs(data[i].M);
    }

    return maxM;
}

double HysteresisAnalyzer::computeRemanence(const std::vector<DataPoint>& data) {
    double minDiff = std::numeric_limits<double>::max();
    double rem = 0.0;

    for(size_t i = 0; i < data.size(); i++) {
        double diff = fabs(data[i].H);
        if(diff < minDiff) {
            minDiff = diff;
            rem = data[i].M;
        }
    }
    return rem;
}

double HysteresisAnalyzer::computeCoercivity(const std::vector<DataPoint>& data) {
    double minDiff = std::numeric_limits<double>::max();
    double Hc = 0.0;

    for(size_t i = 0; i < data.size(); i++) {
        double diff = fabs(data[i].M);
        if(diff < minDiff) {
            minDiff = diff;
            Hc = data[i].H;
        }
    }
    return Hc; 
}

double HysteresisAnalyzer::computeLoopArea(const std::vector<DataPoint>& data) {
    double area = 0.0;
    for (size_t i = 1; i < data.size(); i++) {
        double dH = data[i].H - data[i-1].H;
        double avgM = 0.5 * (data[i].M + data[i-1].M);
        area += avgM * dH;
    }

    return fabs(area);
}