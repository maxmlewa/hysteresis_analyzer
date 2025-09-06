#include "DataLoader.h"
#include "HysteresisAnalyzer.h"
#include "ModelFitter.h"
#include <iostream>

int main(int argc, char* argv[]){
    if (argc < 2) {
        std::cerr << "Usage: hysteresis <datafile> [Ms]" << std::endl;
        return -1;
    }

    std::string filename = argv[1];
    double Ms = 1.0;
    if(argc >= 3) Ms = atof(argv[2]);

    std::vector<DataPoint> data = DataLoader::loadData(filename, Ms);
    if (data.empty()){
        std::cerr << "Error: no data loaded." <<std::endl;
        return 1;
    }

    double Ms_est = HysteresisAnalyzer::computeSaturation(data);
    double Mr = HysteresisAnalyzer::computeRemanence(data);
    double Hc = HysteresisAnalyzer::computeCoercivity(data);
    double area = HysteresisAnalyzer::computeLoopArea(data);

    std::cout << "--- Hysteresis Loop Parameters ---" << std::endl;
    std::cout << "Saturation (Ms est): " << Ms_est << std::endl;
    std::cout << "Remanence (Mr): " << Mr << std::endl;
    std::cout << "Coercivity (Hc): " << Hc << std::endl;
    std::cout << "Loop Area: " << area << std::endl;

    std::array<double, 4> params = ModelFitter::fitModel(data, Ms);
    std::cout << "\n--- Parametric Model Fit ---" << std::endl;
    std::cout << "a = " << params[0] << std::endl;
    std::cout << "b = " << params[1] << std::endl;
    std::cout << "c = " << params[2] << std::endl;
    std::cout << "d = " << params[3] << std::endl;

    return 0;
}