#ifndef PILON_NORMAL_H
#define PILON_NORMAL_H

#include <vector>
#include <string>
#include <cmath>

namespace pilon {

class NormalDistribution {
public:
    // Construct from int array (most common use: coverage, badCoverage, etc.)
    NormalDistribution(const std::vector<int>& values, int nMoments = 2);

    double mean;
    double median;
    std::vector<double> moments;

    double toSigma(double value) const;
    int toSigma10x(double value) const;
    double fromSigma(double sigma) const;

    std::string toString() const;
};

} // namespace pilon

#endif
