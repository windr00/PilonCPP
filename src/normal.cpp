#include "normal.h"
#include <algorithm>
#include <sstream>
#include <cmath>

namespace pilon {

NormalDistribution::NormalDistribution(const std::vector<int>& values, int nMoments)
    : mean(0.0), median(0.0) {
    if (values.empty() || nMoments < 2) return;

    // Compute mean
    double total = 0;
    for (int v : values) total += v;
    mean = total / values.size();

    // Compute median
    std::vector<int> sorted(values);
    std::sort(sorted.begin(), sorted.end());
    int n = sorted.size();
    if (n % 2 == 0)
        median = (sorted[n / 2] + sorted[n / 2 - 1]) / 2.0;
    else
        median = sorted[n / 2];

    // Compute moments (sigma, etc.)
    moments.resize(nMoments);
    moments[0] = mean;
    for (int m = 1; m < nMoments; m++) {
        double sum = 0;
        for (int v : values) {
            sum += std::pow(v - mean, m + 1);
        }
        moments[m] = std::pow(sum / values.size(), 1.0 / (m + 1));
    }
}

double NormalDistribution::toSigma(double value) const {
    if (moments.size() < 2 || moments[1] == 0) return 0;
    return (value - moments[0]) / moments[1];
}

int NormalDistribution::toSigma10x(double value) const {
    return static_cast<int>(std::round(toSigma(value) * 10.0));
}

double NormalDistribution::fromSigma(double sigma) const {
    return moments[0] + sigma * moments[1];
}

std::string NormalDistribution::toString() const {
    std::ostringstream oss;
    oss << "<moments: n=" << moments.size();
    if (!moments.empty()) {
        oss << ",moments=" << moments[0];
        for (size_t i = 1; i < moments.size(); i++)
            oss << "," << moments[i];
    }
    oss << ">";
    return oss.str();
}

} // namespace pilon
