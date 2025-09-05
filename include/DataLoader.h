#ifndef DATA_LOADER_H
#define DATA_LOADER_H

#include <vector>
#include <string>

struct DataPoint {
    double time;    // 1st column (in seconds)
    double H;       // 2nd column (magnetic field in Oersted)
    double M;       // Normalized magnetic signal
};

class DataLoader {
public:
    static std::vector<DataPoint> loadData(const std::string& filename,
                                           double Ms = 1.0);
};

#endif