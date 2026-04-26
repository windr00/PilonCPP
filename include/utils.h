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

#ifndef PILON_UTILS_H
#define PILON_UTILS_H

#include <string>
#include <vector>
#include <cmath>

namespace pilon {

class Utils {
public:
    // Round and divide: (numerator + denominator/2) / denominator
    static int roundDiv(long long numerator, long long denominator);

    // Percentage: (numerator * 100 + denominator/2) / denominator
    static int pct(long long numerator, long long denominator);

    // Split string by delimiter
    static std::vector<std::string> split(const std::string& s, char delimiter);

    // Trim whitespace
    static std::string trim(const std::string& s);

    // Read entire file into string
    static std::string readFile(const std::string& path);

    // Write string to file
    static bool writeFile(const std::string& path, const std::string& content);

    // Check if file exists
    static bool fileExists(const std::string& path);

    // Create directories recursively
    static bool mkdirs(const std::string& path);
};

} // namespace pilon

#endif // PILON_UTILS_H
