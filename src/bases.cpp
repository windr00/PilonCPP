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

#include "bases.h"
#include <algorithm>

namespace pilon {

int Bases::toIndex(char c) {
    switch (c) {
        case 'A': return A;
        case 'C': return C;
        case 'G': return G;
        case 'T': return T;
        default: return -1;
    }
}

char Bases::toBase(int c) {
    switch (c) {
        case A: return 'A';
        case C: return 'C';
        case G: return 'G';
        case T: return 'T';
        default: return 'N';
    }
}

char Bases::toBase(int, int) {
    return 'x';
}

char Bases::complement(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return c;
    }
}

int Bases::complementIndex(int i) {
    return i ^ 3;
}

std::string Bases::reverseComplement(const std::string& s) {
    std::string result(s.size(), ' ');
    for (size_t i = 0; i < s.size(); ++i) {
        result[s.size() - 1 - i] = complement(s[i]);
    }
    return result;
}

int Bases::bit(int i) {
    return 1 << i;
}

const std::unordered_map<char, int>& Bases::getMapToBits() {
    static const std::unordered_map<char, int> mapToBits = {
        {'A', Abit},
        {'C', Cbit},
        {'G', Gbit},
        {'T', Tbit},
        {'R', Abit | Gbit},
        {'Y', Cbit | Tbit},
        {'S', Gbit | Cbit},
        {'W', Abit | Tbit},
        {'K', Gbit | Tbit},
        {'M', Abit | Cbit},
        {'B', Cbit | Gbit | Tbit},
        {'D', Abit | Gbit | Tbit},
        {'H', Abit | Cbit | Tbit},
        {'V', Abit | Cbit | Gbit},
        {'N', Abit | Cbit | Gbit | Tbit}
    };
    return mapToBits;
}

const std::unordered_map<int, char>& Bases::getMapToIUPAC() {
    static const std::unordered_map<int, char> mapToIUPAC = []() {
        std::unordered_map<int, char> m;
        const auto& bits = getMapToBits();
        for (const auto& p : bits) {
            m[p.second] = p.first;
        }
        return m;
    }();
    return mapToIUPAC;
}

char Bases::toIUPAC(char base1, char base2) {
    int i1 = toIndex(base1);
    int i2 = toIndex(base2);
    if (i1 < 0 || i2 < 0) return 'N';
    int bits = bit(i1) | bit(i2);
    const auto& map = getMapToIUPAC();
    auto it = map.find(bits);
    if (it != map.end()) return it->second;
    return 'N';
}

bool Bases::baseMatch(char iupac, char base) {
    if (iupac == base) return true;
    const auto& map = getMapToBits();
    auto iuIt = map.find(iupac);
    auto baseIt = map.find(base);
    if (iuIt == map.end() || baseIt == map.end()) return false;
    return (iuIt->second & baseIt->second) == baseIt->second;
}

} // namespace pilon
