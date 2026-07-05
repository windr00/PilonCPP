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

#ifndef PILON_SCAFFOLD_H
#define PILON_SCAFFOLD_H

#include <string>
#include <vector>
#include <unordered_map>
#include <functional>
#include <cstdint>
#include "assembler.h"
#include "bamfile.h"

namespace pilon {

class MatePair {
public:
    int scaffold1;
    int coord1;
    int scaffold2;
    int coord2;
    int mq;

    MatePair(const BamRead& r1, const BamRead& r2,
             const std::unordered_map<std::string, int>& refNameToIdx);

    bool sameScaffold() const { return scaffold1 == scaffold2; }

    int64_t longCoord1() const;
    int64_t longCoord2() const;
    int64_t distance() const;
    bool isRunt(int maxInsert = 10000) const;
    bool ambiguousPlacement() const { return mq < 4; }
};

class LinkCluster {
public:
    int nLinks;
    int scaffold1;
    int minCoord1;
    int maxCoord1;
    int scaffold2;
    int minCoord2;
    int maxCoord2;
    int mq;
    std::string seqName1;
    std::string seqName2;
    int size1;
    int size2;

    LinkCluster(const std::vector<MatePair>& mates,
                const std::vector<std::string>& refNames,
                const std::vector<int>& refLengths,
                int sigma);

    bool sameScaffold() const { return scaffold1 == scaffold2; }
    bool sameOrientation() const {
        return (minCoord1 < 0 && minCoord2 < 0) || (minCoord1 > 0 && minCoord2 > 0);
    }
    int spread1() const { return maxCoord1 - minCoord1; }
    int spread2() const { return maxCoord2 - minCoord2; }
    bool valid() const { return spread1() > sigma_ && spread2() > sigma_; }
    bool nearEnd(int coord, int size) const { return (coord < sigma_) || (-coord > size - sigma_); }
    bool nearEnds() const { return nearEnd(minCoord1, size1) && nearEnd(minCoord2, size2); }
    bool circular() const { return valid() && sameScaffold() && nearEnds(); }
    bool scaffoldLink() const { return valid() && nearEnds() && !sameScaffold(); }
    bool rearrangement() const { return valid() && !circular() && !scaffoldLink(); }

    std::string reportCoord1() const;
    std::string reportCoord2() const;

private:
    int sigma_;
};

class EndAlignment {
public:
    std::string readName;
    std::string contigName;
    int contigLength;
    int alignStart;
    int alignEnd;
    int unclippedStart;
    int unclippedEnd;
    bool rc;
    bool cantileveredLeft;
    bool cantileveredRight;

    static constexpr int EndMinOverlap = 500;
    static constexpr int EndMaxLength = 5000;

    EndAlignment(const BamRead& read, const std::string& cName, int cLen);
    bool cantilevered() const { return cantileveredLeft || cantileveredRight; }
    std::string contigEndName() const;
};

class EndAlignmentPair {
public:
    EndAlignment ca1;
    EndAlignment ca2;

    EndAlignmentPair(const EndAlignment& a1, const EndAlignment& a2);
    int impliedLength() const;
};

namespace Scaffold {

    static constexpr int EndMinLinks = 5;

    std::vector<std::vector<MatePair>> findClustersInternal(
        const std::vector<MatePair>& coords, int window, int minCluster,
        std::function<int64_t(const MatePair&)> mpFunc);

    std::vector<LinkCluster> findClusters(
        std::vector<MatePair>& coords,
        const std::vector<std::string>& refNames,
        const std::vector<int>& refLengths,
        int sigma, int minCluster = 25);

    void analyzeStrays(BamFile& bam);
    void analyze(std::vector<BamFile*>& bamFiles);

    std::vector<EndAlignment> findEndAlignments(BamFile& bam);
    std::unordered_map<std::string, int> findCircles(const std::vector<EndAlignment>& eas);
    std::unordered_map<std::string, int> findHgapCircles(std::vector<BamFile*>& bamFiles);
    int estimateLength(const std::vector<int>& estimates);

    void dumpCoords(const std::vector<MatePair>& coords);
}

} // namespace pilon

#endif // PILON_SCAFFOLD_H
