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

#ifndef PILON_BASESUM_H
#define PILON_BASESUM_H

#include <array>
#include <string>
#include <vector>

namespace pilon {

class BaseSum {
public:
    std::array<long long, 4> sums;

    BaseSum();

    void add(int base, int n = 1);
    void remove(int base, int n = 1);
    long long sum() const;

    // Returns {maxIndex, max, totalSum}
    std::tuple<int, long long, long long> maxBase() const;

    // Returns {base, excess, total}
    std::tuple<int, long long, long long> excessBase() const;

    // Returns indices sorted by descending count
    std::vector<int> order() const;

    BaseSum operator/(int divisor) const;

    std::string toStringPct() const;
    std::string toString() const;
};

} // namespace pilon

#endif // PILON_BASESUM_H
