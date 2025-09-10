#include "PolynomialModel.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

// Least-squares polynomial fit
static std::vector<double> polyfit(const std::vector<double>& x,
                                   const std::vector<double>& y,
                                   int degree) {
    int n = degree + 1;
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
    std::vector<double> b(n, 0.0);

    for (size_t i = 0; i < x.size(); i++) {
        std::vector<double> xp(n, 1.0);
        for (int j = 1; j < n; j++) xp[j] = xp[j-1] * x[i];
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) A[j][k] += xp[j] * xp[k];
            b[j] += xp[j] * y[i];
        }
    }

    for (int i = 0; i < n; i++) {
        int piv = i;
        for (int r = i+1; r < n; r++)
            if (fabs(A[r][i]) > fabs(A[piv][i])) piv = r;
        std::swap(A[i], A[piv]); std::swap(b[i], b[piv]);

        double diag = A[i][i];
        if (fabs(diag) < 1e-12) continue;
        for (int c = i; c < n; c++) A[i][c] /= diag;
        b[i] /= diag;

        for (int r = 0; r < n; r++) {
            if (r == i) continue;
            double f = A[r][i];
            for (int c = i; c < n; c++) A[r][c] -= f * A[i][c];
            b[r] -= f * b[i];
        }
    }
    return b;
}

PolynomialFitResult PolynomialModel::fit(const std::vector<DataPoint>& data, int degree) {
    PolynomialFitResult r;
    r.degree = degree;

    std::vector<double> Hup, Mup, Hdown, Mdown;
    for (size_t i = 1; i < data.size(); i++) {
        if (data[i].H > data[i-1].H) {
            Hup.push_back(data[i].H);
            Mup.push_back(data[i].M);
        } else {
            Hdown.push_back(data[i].H);
            Mdown.push_back(data[i].M);
        }
    }

    if (!Hup.empty()) r.coeffs_up = polyfit(Hup, Mup, degree);
    if (!Hdown.empty()) r.coeffs_down = polyfit(Hdown, Mdown, degree);

    r.Hmin = std::min_element(data.begin(), data.end(),
                              [](auto &a, auto &b){ return a.H < b.H; })->H;
    r.Hmax = std::max_element(data.begin(), data.end(),
                              [](auto &a, auto &b){ return a.H < b.H; })->H;
    return r;
}

double PolynomialModel::eval(const std::vector<double>& coeffs, double H) {
    double y = 0.0, powH = 1.0;
    for (size_t i = 0; i < coeffs.size(); i++) {
        y += coeffs[i] * powH;
        powH *= H;
    }
    return y;
}

double PolynomialModel::computeSaturation(const PolynomialFitResult& r) {
    return std::max(fabs(eval(r.coeffs_up, r.Hmax)),
                    fabs(eval(r.coeffs_down, r.Hmax)));
}

double PolynomialModel::computeRemanence(const PolynomialFitResult& r) {
    return 0.5 * (eval(r.coeffs_up, 0.0) + eval(r.coeffs_down, 0.0));
}

double PolynomialModel::computeCoercivity(const PolynomialFitResult& r) {
    double Hc_up = 0.0, Hc_down = 0.0;
    double best_up = 1e9, best_down = 1e9;
    for (double H = r.Hmin; H <= r.Hmax; H += (r.Hmax - r.Hmin) / 2000.0) {
        double Mu = eval(r.coeffs_up, H);
        double Md = eval(r.coeffs_down, H);
        if (fabs(Mu) < best_up) { best_up = fabs(Mu); Hc_up = H; }
        if (fabs(Md) < best_down) { best_down = fabs(Md); Hc_down = H; }
    }
    return 0.5 * (Hc_up + Hc_down);
}

double PolynomialModel::computeLoopArea(const PolynomialFitResult& r, int N) {
    double area = 0.0;
    double dH = (r.Hmax - r.Hmin) / N;
    for (int i = 1; i < N; i++) {
        double H1 = r.Hmin + (i-1) * dH;
        double H2 = r.Hmin + i * dH;
        double Mup1 = eval(r.coeffs_up, H1), Mup2 = eval(r.coeffs_up, H2);
        double Mdown1 = eval(r.coeffs_down, H1), Mdown2 = eval(r.coeffs_down, H2);
        area += 0.5 * ((Mup1 - Mdown1) + (Mup2 - Mdown2)) * (H2 - H1);
    }
    return fabs(area);
}

bool PolynomialModel::write_csv(const std::string& path, const PolynomialFitResult& r, int N) {
    std::ofstream f(path.c_str());
    if (!f.is_open()) return false;
    f << "H,Mup,Mdown\n";
    double dH = (r.Hmax - r.Hmin) / N;
    for (int i = 0; i < N; i++) {
        double H = r.Hmin + i * dH;
        double Mu = eval(r.coeffs_up, H);
        double Md = eval(r.coeffs_down, H);
        f << std::setprecision(10) << H << "," << Mu << "," << Md << "\n";
    }
    return true;
}