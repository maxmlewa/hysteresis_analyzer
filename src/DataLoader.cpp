#include "DataLoader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>

static double median_of_vector(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    size_t n = v.size();
    if (n % 2 == 1) return v[n/2];
    else return 0.5*(v[n/2 - 2] + v[n/2]);
    
}

std::vector<DataPoint> DataLoader::loadData(const std::string& filename,
                                            double Ms,
                                            bool doOutlierRemoval,
                                            double outlierK,
                                            const std::string& scaleMode,
                                            double percentile) {
    std::vector<DataPoint> data;
    std::ifstream infile(filename.c_str());
    
    if(!infile.is_open()) {
        std::cerr<< "Error: cannot open " << filename << std::endl;
        return data;
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        double t, H, V, col4, col5;
        if (!(iss >> t >> H >> V >> col4 >> col5)) continue;

        DataPoint p;
        p.time = t;
        p.H = H;
        p.M = V; // temporary, to be processed
        data.push_back(p);
    }

    if (data.empty()){
        std::cerr << "Error: no valid data loaded from " << filename << std::endl;
        return data;
    }

    // Data preprocessing
    
    // mean and stddev
    double meanV = 0.0;
    for (size_t i = 0; i<data.size(); i++) {
        meanV += data[i].M; // M still contains the raw voltage values
    }

    meanV /= data.size();

    double variance = 0.0;
    for (size_t i = 0; i < data.size(); i++) {
        double diff = data[i].M - meanV;
        variance += diff*diff;
    }

    double stddev = std::sqrt(variance / data.size());

    // building the clean values with optional outlier removal
    std::vector<double> cleanVals;
    cleanVals.reserve(data.size());
    for (size_t i = 0; i < data.size(); i++) {
        double v = data[i].M;
        if (!doOutlierRemoval || std::fabs(data[i].M - meanV) <= outlierK * stddev) {
            cleanVals.push_back(data[i].M);
        }
    }

    if( cleanVals.empty()) {
        std::cerr << "Warning : all data points rejected as outliers. Using the full dataset.\n";
        for (size_t i = 0; i < data.size(); i++) {
            cleanVals.push_back(data[i].M);
        }
    }

    // recomputing the mean from the cleaned data
    double cleanMean = 0.0;
    for (size_t i = 0; i < cleanVals.size(); i++) {
        cleanMean += cleanVals[i]; 
    }

    cleanMean /= cleanVals.size();

    // absolute deviations from the cleanMean

    std::vector<double> deviations;
    deviations.reserve(cleanVals.size());
    for (size_t i = 0; i < cleanVals.size(); i++) {
        deviations.push_back(std::fabs(cleanVals[i] - cleanMean));
    }

    // determine the metric scale
    double scaleVal = 0.0;
    if (scaleMode == "max") {
        scaleVal = *std::max_element(deviations.begin(), deviations.end());
    }

    else if(scaleMode == "rms") {
        double s = 0.0;
        for (size_t i = 0; i < deviations.size(); ++i) {
            s += deviations[i] * deviations [i];
        }
        scaleVal = std::sqrt(s/std::max<size_t>(1, deviations.size()));
    }

    else if(scaleMode == "mad") {
        // MAD: Median Absolute Deviation
        scaleVal = median_of_vector(deviations);

        // for MAD to be comparable to the normal distribution, multiply by 1.4826
        scaleVal *= 1.4826;
    }

    else if(scaleMode == "percentile") {
        std::sort(deviations.begin(), deviations.end());
        double p = percentile;
        if(p <= 0.0) p = 0.01;
        if(p >= 1.0) p = 0.99;

        size_t idx = static_cast<size_t>(std::floor(p * (deviations.size() - 1)));
        scaleVal = deviations[idx];
    }

    else { // RMS fallback
        double s = 0.0;
        for (size_t i = 0; i < deviations.size(); i++) {
            s += deviations[i] * deviations[i];
        }
        scaleVal = std::sqrt(s/ std::max<size_t>(1, deviations.size()));

    }

    if (scaleVal <= 0.0) {
        std::cerr << "Warning: scale metric is zero or tiny, using scale = 1.0\n";
        scaleVal = 1.0;
    }

    double scale = Ms/scaleVal;

    // applying correction to all datapoints
    for (size_t i = 0; i < data.size(); ++i){
        data[i].M = (data[i].M - cleanMean) * scale;
    }

    return data;
                                            }