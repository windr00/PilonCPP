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

#ifndef PILON_PILEUP_H
#define PILON_PILEUP_H

#include <string>
#include <vector>
#include <utility>
#include "basesum.h"

namespace pilon {

// Forward declaration
class Pilon;

class PileUp {
public:
    BaseSum baseCount;
    BaseSum qualSum;
    long long mqSum;
    long long qSum;
    long long physCov;
    long long insertSize;
    long long badPair;
    long long deletions;
    long long delQual;
    long long insertions;
    long long insQual;
    long long clips;
    std::vector<std::vector<uint8_t>> insertionList;
    std::vector<std::vector<uint8_t>> deletionList;

    PileUp();

    long long count() const;
    long long depth() const;

    int baseIndex(char c) const;
    char indexBase(int i) const;

    int weightedMq() const;
    int weightedQual() const;
    int meanQual() const;
    int meanMq() const;

    void add(char base, int qual, int mq);
    void remove(char base, int qual, int mq);

    void addInsertion(const std::vector<uint8_t>& insertion, int qual, int mq);
    void addDeletion(const std::vector<uint8_t>& deletion, int qual, int mq);

    long long totalQSum() const;

    int insPct() const;
    int delPct() const;
    int qualPct() const;
    int clipPct() const;

    class BaseCall {
    public:
        long long n;
        int baseIndex;
        int altBaseIndex;
        char base;
        long long baseSum;
        char altBase;
        long long altBaseSum;
        bool homo;
        int score;
        std::string insertion;
        std::string deletion;
        bool indel;
        bool homoIndel;

        BaseCall(const PileUp& pu);

        bool isInsertion() const;
        bool isDeletion() const;
        bool called() const;
        int q() const;
        bool highConfidence() const;

        std::string callString(bool indelOk = true) const;
        bool baseMatch(char refBase) const;
        char iupacBase() const;

        std::string toString() const;
    private:
        std::string insertCall() const;
        std::string deletionCall() const;
        std::string indelCall(const std::vector<std::vector<uint8_t>>& indelList, int pct) const;
        std::pair<std::string, bool> hetIndelCall(const std::vector<std::vector<uint8_t>>& indelList, int pct) const;
    };

    BaseCall baseCall() const;

    std::string toString() const;
};

} // namespace pilon

#endif // PILON_PILEUP_H
