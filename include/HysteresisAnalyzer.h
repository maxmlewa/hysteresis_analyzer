#ifndef HYSTERESIS_ANALYZER_H
#define HYSTERESIS_ANALYZER_H

#include "DataLoader.h"
#include <vector>

class HysteresisAnalyzer {
public:
    static double computeSaturation(const std::vector<DataPoint>& data);
    static double computeRemanence(const std::vector<DataPoint>& data);
    static double computeCoercivity(const std::vector<DataPoint>& data);
    static double computeLoopArea(const std::vector<DataPoint>& data);
};

#endif