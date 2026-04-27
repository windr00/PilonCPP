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

#ifndef PILON_GAPFILLER_H
#define PILON_GAPFILLER_H

#include <string>
#include <vector>
#include <utility>
#include "assembler.h"
#include "genome.h"

namespace pilon {

class BamFile;

class GapFiller {
public:
    static constexpr int minExtend = 20;
    static int k;  // 2 * Assembler::K + 1

    GapFiller(const GenomeRegion& region);

    // Returns {start, refSeq, patchSeq} or {0, "", ""} for no solution
    std::tuple<int, std::string, std::string> fillGap(const Region& gap);
    std::tuple<int, std::string, std::string> fixBreak(const Region& brk);

    std::string tandemRepeat() const { return tandemRepeat_; }

    // Static methods called from GenomeRegion
    static void fixBreak(GenomeRegion& region, const Region& breakRegion);
    static void fixGap(GenomeRegion& region, const Region& gapRegion);
    static bool fixCircles(GenomeRegion& region, int gapMargin);
    static void fixNovel(GenomeFile* genome, const std::vector<BamFile*>& bamFiles);

private:
    const GenomeRegion& region_;
    std::string tandemRepeat_;
    static const std::tuple<int, std::string, std::string> noSolution;

    std::tuple<int, std::string, std::string>
    assembleAcrossBreak(const Region& brk, bool isGap);

    std::string consensusFromLeft(const std::vector<std::string>& seqs);
    std::string consensusFromRight(const std::vector<std::string>& seqs);

    std::tuple<int, std::vector<std::string>, std::vector<std::string>, int, std::string>
    assembleIntoBreak(const Region& brk, const std::vector<BamRead>& reads);

    std::vector<std::tuple<int, std::string, std::string>>
    breakJoins(int start, const std::vector<std::string>& forwardPaths,
               const std::vector<std::string>& reversePaths, int stop);

    std::tuple<int, std::string, std::string>
    joinBreak(int startArg, const std::string& forward,
              const std::string& reverse, int stopArg);

    std::string properOverlap(const std::string& left, const std::string& right,
                              int minOverlap);

    std::tuple<int, std::string, std::string>
    trimPatch(int startArg, std::string patchArg, int stopArg);

    bool partialMatchesReference(int start, const std::string& fromLeft,
                                 const std::string& fromRight, int stop,
                                 int loopLength);

    int breakRadius() const;

    std::vector<BamRead> recruitReads(const Region& brk);
};

} // namespace pilon

#endif // PILON_GAPFILLER_H
