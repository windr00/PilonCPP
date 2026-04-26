/*
 * Copyright (c) 2012-2018 Broad Institute, Inc.
 *
 * This file is part of PilonCpp.
 *
 * PilonCpp is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2
 * as published by the Free Software Foundation.
 *
 * PilonCpp is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PilonCpp.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "basesum.h"
#include <numeric>
#include <algorithm>
#include <sstream>

namespace pilon {

BaseSum::BaseSum() : sums{0, 0, 0, 0} {}

void BaseSum::add(int base, int n) {
    if (base >= 0 && base < 4) sums[base] += n;
}

void BaseSum::remove(int base, int n) {
    if (base >= 0 && base < 4) sums[base] -= n;
}

long long BaseSum::sum() const {
    return std::accumulate(sums.begin(), sums.end(), 0LL);
}

std::tuple<int, long long, long long> BaseSum::maxBase() const {
    int maxIndex = 0;
    long long max = sums[0];
    long long total = max;
    for (int b = 1; b < 4; b++) {
        total += sums[b];
        if (sums[b] > max) {
            max = sums[b];
            maxIndex = b;
        }
    }
    return {maxIndex, max, total};
}

std::tuple<int, long long, long long> BaseSum::excessBase() const {
    auto [base, max, total] = maxBase();
    long long excess = max;
    for (int b = 0; b < 4; b++) {
        if (b != base) excess -= sums[b];
    }
    return {base, excess, total};
}

std::vector<int> BaseSum::order() const {
    std::vector<int> indices = {0, 1, 2, 3};
    std::sort(indices.begin(), indices.end(),
              [this](int a, int b) { return sums[a] > sums[b]; });
    return indices;
}

BaseSum BaseSum::operator/(int divisor) const {
    BaseSum result;
    if (divisor != 0) {
        for (int i = 0; i < 4; i++) {
            result.sums[i] = sums[i] / divisor;
        }
    }
    return result;
}

std::string BaseSum::toStringPct() const {
    long long div = sum();
    std::ostringstream oss;
    for (int i = 0; i < 4; i++) {
        if (i > 0) oss << ",";
        oss << (div == 0 ? 0 : (100 * sums[i] + div / 2) / div);
    }
    return oss.str();
}

std::string BaseSum::toString() const {
    std::ostringstream oss;
    for (int i = 0; i < 4; i++) {
        if (i > 0) oss << ",";
        oss << sums[i];
    }
    return oss.str();
}

} // namespace pilon
