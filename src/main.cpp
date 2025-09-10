#include "DataLoader.h"
#include "HysteresisAnalyzer.h"
#include "ModelFitter.h"
#include "ParametricModel.h"
#include "PolynomialModel.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <filesystem>

namespace fs = std::filesystem;

// ----------------- CSV Writers -----------------

// Preprocessed data CSV
static bool write_data_csv(const std::string &path, const std::vector<DataPoint> &data) {
    std::ofstream f(path.c_str());
    if (!f.is_open()) return false;
    f << "H,M\n";
    for (const auto &d : data) f << std::setprecision(10) << d.H << "," << d.M << "\n";
    return true;
}

// Parametric model CSV
static bool write_model_csv(const std::string &path, const ParametricModel::Params &p, size_t N=2000) {
    if (p.empty()) return false;
    std::vector<double> Hs, Ms;
    ParametricModel::sample_curve(p, N, Hs, Ms);
    std::ofstream f(path.c_str());
    if (!f.is_open()) return false;
    f << "H,ModelM\n";
    for (size_t i = 0; i < Hs.size(); ++i) f << std::setprecision(10) << Hs[i] << "," << Ms[i] << "\n";
    return true;
}

// Write metrics JSON (simple)
static bool write_metrics_json(const std::string &path,
                               const FitResult &fit,
                               double dataMs, double dataMr, double dataHc, double dataArea,
                               double polyMs, double polyMr, double polyHc, double polyArea) {
    std::ofstream f(path.c_str());
    if (!f.is_open()) return false;
    f << "{\n";
    f << "  \"data\": {\n";
    f << "    \"Ms\": " << dataMs << ",\n";
    f << "    \"Mr\": " << dataMr << ",\n";
    f << "    \"Hc\": " << dataHc << ",\n";
    f << "    \"Area\": " << dataArea << "\n";
    f << "  },\n";
    f << "  \"parametric\": {\n";
    f << "    \"Ms\": " << fit.Ms_model << ",\n";
    f << "    \"Mr\": " << fit.Mr_model << ",\n";
    f << "    \"Hc\": " << fit.Hc_model << ",\n";
    f << "    \"Area\": " << fit.area_model << "\n";
    f << "  },\n";
    f << "  \"polynomial\": {\n";
    f << "    \"Ms\": " << polyMs << ",\n";
    f << "    \"Mr\": " << polyMr << ",\n";
    f << "    \"Hc\": " << polyHc << ",\n";
    f << "    \"Area\": " << polyArea << "\n";
    f << "  }\n";
    f << "}\n";
    return true;
}

static void print_usage_and_exit() {
    std::cout <<
    "Usage: hysteresis <datafile> [options]\n"
    "Options:\n"
    "  --Ms VAL              target Ms normalization (default 1.0)\n"
    "  --plot [none|data|model|both]\n"
    "  --plot-auto           run Python plot script automatically (requires python3+matplotlib)\n"
    "  --outdir PATH         output folder for CSV/JSON (default ./output)\n";
    exit(1);
}

int main(int argc, char* argv[]) {
    if (argc < 2) print_usage_and_exit();

    std::string filename = argv[1];
    double Ms = 1.0;
    std::string plotMode = "both";
    bool plotAuto = false;
    std::string outdir = "output";

    // parse args
    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--Ms" && i+1 < argc) { Ms = atof(argv[++i]); }
        else if (a == "--plot" && i+1 < argc) { plotMode = argv[++i]; }
        else if (a == "--plot-auto") { plotAuto = true; }
        else if (a == "--outdir" && i+1 < argc) { outdir = argv[++i]; }
        else { std::cerr << "Unknown arg: " << a << "\n"; print_usage_and_exit(); }
    }

    try {
        fs::create_directories(outdir);
    } catch (const std::exception& e) {
        std::cerr << "Error creating output directory: " << e.what() << "\n";
        return 1;
    }

    // Load & preprocess data
    auto data = DataLoader::loadData(filename, Ms, true, 3.0, "rms", 0.95);
    if (data.empty()) { std::cerr << "Failed to load data\n"; return 1; }

    // Analyzer metrics
    double dataMs = HysteresisAnalyzer::computeSaturation(data);
    double dataMr = HysteresisAnalyzer::computeRemanence(data);
    double dataHc = HysteresisAnalyzer::computeCoercivity(data);
    double dataArea = HysteresisAnalyzer::computeLoopArea(data);

    // Parametric fit
    double maxH = 0.0;
    for (auto &d : data) maxH = std::max(maxH, std::fabs(d.H));
    ParametricModel::Params guess(10, 0.0);
    guess[0] = 0.9; guess[1] = 1.0; guess[2] = 1.0;
    guess[3] = 3.0; guess[4] = 1.0;
    guess[5] = guess[6] = guess[7] = 0.0;
    guess[8] = maxH; guess[9] = dataMs;

    FitResult res;
    if (!ModelFitter::fit_staged(data, guess, res)) {
        std::cerr << "Parametric model fit failed\n";
    }

    // Polynomial fit
    PolynomialFitResult polyRes = PolynomialModel::fit(data, 9);
    double polyMs = PolynomialModel::computeSaturation(polyRes);
    double polyMr = PolynomialModel::computeRemanence(polyRes);
    double polyHc = PolynomialModel::computeCoercivity(polyRes);
    double polyArea = PolynomialModel::computeLoopArea(polyRes);

    // Tabular output
    std::cout << "\n=== Hysteresis Loop Metrics ===\n";
    std::cout << std::setw(12) << "Source"
              << std::setw(12) << "Ms"
              << std::setw(12) << "Mr"
              << std::setw(12) << "Hc"
              << std::setw(12) << "Area" << "\n";
    std::cout << std::string(60, '-') << "\n";

    std::cout << std::setw(12) << "Data"
              << std::setw(12) << dataMs
              << std::setw(12) << dataMr
              << std::setw(12) << dataHc
              << std::setw(12) << dataArea << "\n";

    std::cout << std::setw(12) << "Parametric"
              << std::setw(12) << res.Ms_model
              << std::setw(12) << res.Mr_model
              << std::setw(12) << res.Hc_model
              << std::setw(12) << res.area_model << "\n";

    std::cout << std::setw(12) << "Polynomial"
              << std::setw(12) << polyMs
              << std::setw(12) << polyMr
              << std::setw(12) << polyHc
              << std::setw(12) << polyArea << "\n";

    // Exports
    std::string data_csv = outdir + "/data_preprocessed.csv";
    std::string model_csv = outdir + "/model_fit.csv";
    std::string poly_csv = outdir + "/poly_fit.csv";
    std::string metrics_json = outdir + "/metrics.json";

    write_data_csv(data_csv, data);
    write_model_csv(model_csv, res.final_params, 3000);
    PolynomialModel::write_csv(poly_csv, polyRes, 2000);
    write_metrics_json(metrics_json, res, dataMs, dataMr, dataHc, dataArea,
                       polyMs, polyMr, polyHc, polyArea);

    // Plotting
    if (plotAuto) {
        std::ostringstream pcmd;
        pcmd << "python3 scripts/plot_hysteresis.py"
             << " --data " << data_csv
             << " --model " << model_csv
             << " --poly " << poly_csv
             << " --metrics " << metrics_json
             << " --which " << plotMode;
        std::cout << "Running: " << pcmd.str() << "\n";
        std::system(pcmd.str().c_str());
    }

    return 0;
}