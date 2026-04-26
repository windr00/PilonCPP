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

#ifndef PILON_BASES_H
#define PILON_BASES_H

#include <string>
#include <unordered_map>

namespace pilon {

class Bases {
public:
    // values chosen such that x^3 is the complementary base
    static constexpr int A = 0;
    static constexpr int C = 1;
    static constexpr int G = 2;
    static constexpr int T = 3;

    static int toIndex(char c);
    static char toBase(int c);
    static char toBase(int c1, int c2);
    static char complement(char c);
    static int complementIndex(int i);
    static std::string reverseComplement(const std::string& s);

    static int bit(int i);
    static constexpr int Abit = 1 << A;
    static constexpr int Cbit = 1 << C;
    static constexpr int Gbit = 1 << G;
    static constexpr int Tbit = 1 << T;

    static const std::unordered_map<char, int>& getMapToBits();
    static const std::unordered_map<int, char>& getMapToIUPAC();

    static char toIUPAC(char base1, char base2);
    static bool baseMatch(char iupac, char base);
};

} // namespace pilon

#endif // PILON_BASES_H
