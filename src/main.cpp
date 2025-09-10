#include "DataLoader.h"
#include "HysteresisAnalyzer.h"
#include "ModelFitter.h"
#include "ParametricModel.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <cstdio>


// Write preprocessed data CSV
static bool write_data_csv(const std::string &path, const std::vector<DataPoint> &data) {
    std::ofstream f(path.c_str());
    if (!f.is_open()) return false;
    f << "H,M\n";
    for (const auto &d : data) f << std::setprecision(10) << d.H << "," << d.M << "\n";
    f.close();
    return true;
}


// Write sampled model CSV (H,ModelM)
static bool write_model_csv(const std::string &path, const ParametricModel::Params &p, size_t N=2000) {
    std::vector<double> Hs, Ms;
    ParametricModel::sample_curve(p, N, Hs, Ms);
    std::ofstream f(path.c_str());
    if (!f.is_open()) return false;
    f << "H,ModelM\n";
    for (size_t i = 0; i < Hs.size(); ++i) f << std::setprecision(10) << Hs[i] << "," << Ms[i] << "\n";
    f.close();
    return true;
}

// Write metrics JSON (no dependency)
static bool write_metrics_json(const std::string &path, const FitResult &fit,
                               double dataMs, double dataMr, double dataHc, double dataArea) {
    std::ofstream f(path.c_str());
    if (!f.is_open()) return false;
    f << "{\n";
    f << "  \"data\": {\n";
    f << "    \"Ms\": " << dataMs << ",\n";
    f << "    \"Mr\": " << dataMr << ",\n";
    f << "    \"Hc\": " << dataHc << ",\n";
    f << "    \"Area\": " << dataArea << "\n";
    f << "  },\n";
    f << "  \"model\": {\n";
    f << "    \"Ms\": " << fit.Ms_model << ",\n";
    f << "    \"Mr\": " << fit.Mr_model << ",\n";
    f << "    \"Hc\": " << fit.Hc_model << ",\n";
    f << "    \"Area\": " << fit.area_model << "\n";
    f << "  }\n";
    f << "}\n";
    f.close();
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

    // parse simple args
    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--Ms" && i+1 < argc) { Ms = atof(argv[++i]); }
        else if (a == "--plot" && i+1 < argc) { plotMode = argv[++i]; }
        else if (a == "--plot-auto") { plotAuto = true; }
        else if (a == "--outdir" && i+1 < argc) { outdir = argv[++i]; }
        else { std::cerr << "Unknown arg: " << a << "\n"; print_usage_and_exit(); }
    }

    // create output dir
    std::string cmd = "mkdir -p " + outdir;
    int rc = std::system(cmd.c_str());
    (void)rc;

    // Load & preprocess data
    auto data = DataLoader::loadData(filename, Ms, true, 3.0, "rms", 0.95);
    if (data.empty()) { std::cerr << "Failed to load data\n"; return 1; }

    // Simple analyzer
    double dataMs = HysteresisAnalyzer::computeSaturation(data);
    double dataMr = HysteresisAnalyzer::computeRemanence(data);
    double dataHc = HysteresisAnalyzer::computeCoercivity(data);
    double dataArea = HysteresisAnalyzer::computeLoopArea(data);

    std::cout << "Data-derived: Ms=" << dataMs << " Mr=" << dataMr
              << " Hc=" << dataHc << " Area=" << dataArea << "\n";

    // Initial guess
    ParametricModel::Params guess = {0.0, 1.0, 1.0, 3.0, 1.0,
                                     0.0, 0.0, 0.0, 1.0, 1.0};
    FitResult res;
    if (!ModelFitter::fit_staged(data, guess, res)) {
        std::cerr << "Model fit failed\n";
    } else {
        std::cout << "Model fit SSE=" << res.sse << "\n";
        std::cout << "Model-derived: Ms=" << res.Ms_model
                  << " Mr=" << res.Mr_model
                  << " Hc=" << res.Hc_model
                  << " Area=" << res.area_model << "\n";
    }

    // Exports
    std::string data_csv = outdir + "/data_preprocessed.csv";
    std::string model_csv = outdir + "/model_fit.csv";
    std::string metrics_json = outdir + "/metrics.json";

    if (!write_data_csv(data_csv, data)) {
        std::cerr << "Warning: failed to write " << data_csv << "\n";
    }
    if (!write_model_csv(model_csv, res.final_params, 3000)) {
        std::cerr << "Warning: failed to write " << model_csv << "\n";
    }
    if (!write_metrics_json(metrics_json, res, dataMs, dataMr, dataHc, dataArea)) {
        std::cerr << "Warning: failed to write " << metrics_json << "\n";
    }

    // Optionally call the Python plotter automatically
    if (plotAuto) {
        std::ostringstream pcmd;
        pcmd << "python3 scripts/plot_hysteresis.py"
             << " --data " << data_csv
             << " --model " << model_csv
             << " --metrics " << metrics_json
             << " --which " << plotMode;
        std::cout << "Running: " << pcmd.str() << "\n";
        int r = std::system(pcmd.str().c_str());
        if (r != 0) {
            std::cerr << "Warning: plotting command returned " << r << "\n";
        }
    }

    std::cout << "Exports written to: " << outdir << "\n";
    std::cout << "To visualize, run:\n"
              << "  python3 scripts/plot_hysteresis.py --data " << data_csv
              << " --model " << model_csv << " --metrics " << metrics_json
              << " --which " << plotMode << "\n";

    return 0;
}