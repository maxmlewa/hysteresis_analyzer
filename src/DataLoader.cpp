#include "DataLoader.h"
#include <fstream>
#include <sstream>
#include <iostream>

std::vector<DataPoint> DataLoader::loadData(const std::string& filename,
                                            double Ms) {
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
        p.M = V/Ms;
        data.push_back(p);
    }

    return data;
                                            }