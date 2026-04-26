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

#include "gapfiller.h"
#include "bases.h"
#include "pilon.h"
#include "bamfile.h"
#include <algorithm>
#include <iostream>
#include <sstream>

namespace pilon {

int GapFiller::k = 2 * Assembler::K + 1;
const std::tuple<int, std::string, std::string> GapFiller::noSolution = {0, "", ""};

GapFiller::GapFiller(const GenomeRegion& region)
    : region_(region) {}

std::tuple<int, std::string, std::string> GapFiller::fillGap(const Region& gap) {
    return assembleAcrossBreak(gap, true);
}

std::tuple<int, std::string, std::string> GapFiller::fixBreak(const Region& brk) {
    return assembleAcrossBreak(brk, false);
}

std::tuple<int, std::string, std::string>
GapFiller::assembleAcrossBreak(const Region& brk, bool isGap) {
    std::vector<BamRead> reads = recruitReads(brk);
    if (reads.empty()) return noSolution;

    auto result = assembleIntoBreak(brk, reads);
    int start = std::get<0>(result);
    auto forwardPaths = std::get<1>(result);
    auto reversePaths = std::get<2>(result);
    int stop = std::get<3>(result);
    std::string loopSeq = std::get<4>(result);

    if (forwardPaths.empty() || reversePaths.empty()) {
        return noSolution;
    }

    // Get consensus sequences
    std::string fromLeft = consensusFromLeft(forwardPaths);
    std::string fromRight = consensusFromRight(reversePaths);

    if (fromLeft.empty() || fromRight.empty()) {
        return noSolution;
    }

    // Check for tandem repeats
    if (!loopSeq.empty()) {
        tandemRepeat_ = loopSeq;
        if (Pilon::trSafe) return noSolution;
    }

    // Find proper overlap
    int minOverlap = isGap ? 10 : 5;
    std::string overlap = properOverlap(fromLeft, fromRight, minOverlap);

    if (overlap.empty()) {
        return noSolution;
    }

    // Join the break
    auto joinResult = joinBreak(start, fromLeft, fromRight, stop);
    return trimPatch(start, std::get<2>(joinResult), stop);
}

std::string GapFiller::consensusFromLeft(const std::vector<std::string>& seqs) {
    if (seqs.empty()) return "";

    // Simple consensus: take the longest sequence
    std::string best = seqs[0];
    for (const auto& s : seqs) {
        if (s.length() > best.length()) {
            best = s;
        }
    }
    return best;
}

std::string GapFiller::consensusFromRight(const std::vector<std::string>& seqs) {
    if (seqs.empty()) return "";

    std::string best = seqs[0];
    for (const auto& s : seqs) {
        if (s.length() > best.length()) {
            best = s;
        }
    }
    return best;
}

std::tuple<int, std::vector<std::string>, std::vector<std::string>, int, std::string>
GapFiller::assembleIntoBreak(const Region& brk, const std::vector<BamRead>& reads) {
    Assembler assembler;
    assembler.addReads(reads);
    assembler.buildGraph();

    // Get sequences from left and right of break
    std::string leftSeq = region_.subString(brk.start - k, k);
    std::string rightSeq = region_.subString(brk.stop, k);

    auto forwardPaths = assembler.tryForward(leftSeq);
    auto reversePaths = assembler.tryReverse(rightSeq);

    int start = brk.start;
    int stop = brk.stop;
    std::string loopSeq = assembler.toString(); // Simplified

    return {start, forwardPaths, reversePaths, stop, loopSeq};
}

std::vector<std::tuple<int, std::string, std::string>>
GapFiller::breakJoins(int start, const std::vector<std::string>& forwardPaths,
                      const std::vector<std::string>& reversePaths, int stop) {
    std::vector<std::tuple<int, std::string, std::string>> joins;

    for (const auto& fwd : forwardPaths) {
        for (const auto& rev : reversePaths) {
            auto result = joinBreak(start, fwd, rev, stop);
            if (!std::get<2>(result).empty()) {
                joins.push_back(result);
            }
        }
    }

    return joins;
}

std::tuple<int, std::string, std::string>
GapFiller::joinBreak(int startArg, const std::string& forward,
                     const std::string& reverse, int stopArg) {
    // Find overlap between forward and reverse
    int minOverlap = 10;
    std::string overlap = properOverlap(forward, reverse, minOverlap);

    if (overlap.empty()) {
        return {0, "", ""};
    }

    // Create patch
    std::string patch = forward + reverse.substr(overlap.length());
    return {startArg, overlap, patch};
}

std::string GapFiller::properOverlap(const std::string& left, const std::string& right,
                                     int minOverlap) {
    int maxOverlap = std::min(static_cast<int>(left.length()),
                              static_cast<int>(right.length()));

    for (int overlap = maxOverlap; overlap >= minOverlap; overlap--) {
        std::string leftEnd = left.substr(left.length() - overlap);
        std::string rightStart = right.substr(0, overlap);

        if (leftEnd == rightStart) {
            return leftEnd;
        }
    }

    return "";
}

std::tuple<int, std::string, std::string>
GapFiller::trimPatch(int startArg, std::string patchArg, int stopArg) {
    if (patchArg.empty()) return {0, "", ""};

    // Trim patch to fit within region
    int maxLen = stopArg - startArg;
    if (static_cast<int>(patchArg.length()) > maxLen) {
        patchArg = patchArg.substr(0, maxLen);
    }

    return {startArg, "", patchArg};
}

bool GapFiller::partialMatchesReference(int start, const std::string& fromLeft,
                                        const std::string& fromRight, int stop,
                                        int loopLength) {
    // Check if the assembled sequence matches the reference
    std::string refSeq = region_.subString(start, stop - start);

    if (fromLeft.length() > refSeq.length()) return false;
    if (fromRight.length() > refSeq.length()) return false;

    // Simple check: see if fromLeft matches beginning of ref
    if (fromLeft.length() > 0) {
        std::string refStart = refSeq.substr(0, fromLeft.length());
        if (fromLeft != refStart) return false;
    }

    return true;
}

int GapFiller::breakRadius() const {
    // Return the radius for searching around breaks
    return 2 * Assembler::K + 100;
}

std::vector<BamRead> GapFiller::recruitReads(const Region& brk) {
    // Recruit reads that span the break/gap region from all BAM files
    std::vector<BamRead> reads;
    
    // Expand the region to capture reads that span the break
    int flank = Pilon::flank;
    Region expanded(brk.name, 
                    std::max(brk.start - flank, 0), 
                    brk.stop + flank);
    
    // Query each BAM file for reads in the expanded region
    for (auto* bam : Pilon::bamFiles) {
        if (!bam) continue;
        
        // Get reads in the flank region
        auto flankReads = bam->recruitFlankReads(expanded);
        
        // Filter to reads that actually span or overlap the break
        for (const auto& read : flankReads) {
            if (read.unmapped || read.secondary || read.duplicate) continue;
            if (read.failsVendorQC && !Pilon::nonPf) continue;
            
            // Check if read overlaps the break region
            bool overlaps = (read.alignmentStart < brk.stop) && 
                           (read.alignmentEnd > brk.start);
            
            if (overlaps) {
                reads.push_back(read);
            }
        }
    }
    
    if (Pilon::verbose && !reads.empty()) {
        std::cout << "  Recruited " << reads.size() 
                  << " reads spanning " << brk.toString() << std::endl;
    }
    
    return reads;
}


} // namespace pilon
