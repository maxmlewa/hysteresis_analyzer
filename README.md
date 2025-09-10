# Parametric Hysteresis Model Fitter (Lapshin Model)

This project provides a C++11 implementation of the exact Lapshin hysteresis model and a Levenberg–Marquardt-based fitting framework for analyzing magnetic hysteresis loops. It extracts parameters, computes residuals, and derives key magnetic properties like saturation magnetization (Ms), remanence (Mr), coercivity (Hc), and loop area. Python scripts are included for visualizing fits. The framework is designed to allow expansion to other models like Jiles–Atherton and Stoner–Wohlfarth.

## Features

- **Exact Lapshin Model**: Implements equations from the Lapshin paper for precise hysteresis loop modeling.
- **Levenberg–Marquardt Optimizer**: Plain C++11 implementation for fitting parameters to experimental data.
- **Model Fitter**: Automatically fits parameters `a, bx, by, m, n, Delta_alpha1, Delta_alpha2, Delta_alpha3` and optional scaling parameters `Hscale` and `Mscale`. Computes residuals and derived magnetic properties.
- **Python Plotting Scripts**: Visualizes the fitted model against experimental data, saving plots in the `output/` folder.
- **Expandable Framework**: Easily integrate other hysteresis models such as Jiles–Atherton and Stoner–Wohlfarth for comparative analysis.

## Build Instructions

1. Clone the repository using `git clone <repo-url>` and navigate to the project folder using `cd <repo-directory>`.
2. Create a build directory with `mkdir -p build` and navigate into it: `cd build`.
3. Run `cmake ..` to configure and `cmake --build .` to build the project.
4. The executable will be located at `build/Debug/hysteresis_analyzer.exe`.

## Quickstart

From the project root, run:

`.\build\Debug\hysteresis_analyzer.exe data/easy_45deg.txt --Ms 1.0 --plot-auto --plot both`

Example output:

`Data-derived: Ms=1.22875 Mr=0.61544 Hc=1113.62 Area=3194.45`
`[ModelFitter] Stage 1: fit scale + basic shape on 300 pts`
`Stage 1 SSE = 34.6112 (ok=true)`
`[ModelFitter] Stage 2: fit exponents (m,n)`
`Stage 2 SSE = 34.6112 (ok=true)`
`[ModelFitter] Stage 3: fit phase shifts`
`Stage 3 SSE = 34.6112 (ok=true)`
`Final SSE = 34.6112 (ok=true)`
`Model fit SSE=34.6112`
`Model-derived: Ms=5.47078 Mr=0.173273 Hc=1229.02 Area=15842.3`
`Running: python3 scripts/plot_hysteresis.py ...`
`Exports written to: output`

## Plotting Results

Run the plotting script manually from the project root:

`python3 scripts/plot_hysteresis.py --data output/data_preprocessed.csv --model output/model_fit.csv --metrics output/metrics.json --which both`

Plots will be saved in the `output/` folder and compare experimental data with the fitted model.

## Notes & Recommendations

- **Identifiability**: Some parameters may be correlated; use physical intuition for initial guesses.
- **Exponents**: Integer exponents `m` and `n` improve numerical stability.
- **Scaling**: Optional `Hscale` and `Mscale` parameters help when data uses non-standard units.
- **Dependencies**: Core fitting framework is pure C++11; Python 3 with matplotlib is only needed for plotting.

## Future Extensions

To integrate new hysteresis models:

1. Create a new C++ class in `include/` and `src/` with methods to set parameters and compute magnetization M(H).
2. Add a fitting routine in `ModelFitter` using the new model class.
3. Use the existing Python scripts to plot and compare results with the Lapshin model.

This allows multi-model comparisons with models like Jiles–Atherton and Stoner–Wohlfarth.

## References

- Lapshin, V. I., *Exact Hysteresis Loop Model for Magnetic Materials*, eqs. (13), (16), (19), (21), (23)
- Jiles–Atherton model: widely used for ferromagnetic materials
- Stoner–Wohlfarth model: exact analytical formulation for hysteresis loops
