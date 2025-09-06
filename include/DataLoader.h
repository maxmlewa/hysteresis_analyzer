#ifndef DATA_LOADER_H
#define DATA_LOADER_H

#include <vector>
#include <string>

struct DataPoint {
    double time;    // 1st column (in seconds)
    double H;       // 2nd column (magnetic field in Oersted)
    double M;       // Normalized magnetic signal
};

// Responsible for loading and preprocessing hysteresis loop data
class DataLoader {
public:
    // Load and normalize hysteresis data.
    // Parameters:          
    //  Ms                  -> target saturation value after scaling (default 1.0)
    //  doOutlierRemoval    -> enable/disable outlier rejection (default true)
    //  outlierK            -> threshold multiplier for sigma clipping
    //  scaleMode           -> "max", "rms", "mad", "percentile" 
    //  percentile          -> used when scaleMode = "percentile"
    static std::vector<DataPoint> loadData( const std::string& filename,
                                            double Ms = 1.0,
                                            bool doOutlierRemoval = true,
                                            double outlierK = 3.0,
                                            const std::string& scaleMode = "rms",
                                            double percentile = 0.95);
};

#endif