/*
 * Copyright (c) 2012-2018 Broad Institute, Inc.
 *
 * This file is part of PilonCpp.
 */

#include "gapfiller.h"
#include "genome.h"
#include "pilon.h"
#include "bamfile.h"
#include "assembler.h"
#include <algorithm>
#include <iostream>
#include <set>
#include <map>
#include <cmath>
#include <tuple>

namespace pilon {

int GapFiller::k = 2 * Assembler::K + 1;
const std::tuple<int, std::string, std::string> GapFiller::noSolution = {0, "", ""};

GapFiller::GapFiller(GenomeRegion& region, std::vector<BamFile*>* bamHandles)
    : region_(region), threadBams_(bamHandles) {}

// Use thread-local BAM handles if available, otherwise fall back to global
#define BAMS() (threadBams_ ? *threadBams_ : Pilon::bamFiles)

static bool substrEq(const std::string& a, int aOff, const std::string& b, int bOff, int len) {
    if (aOff + len > (int)a.length() || bOff + len > (int)b.length()) return false;
    for (int i = 0; i < len; i++)
        if (a[aOff + i] != b[bOff + i]) return false;
    return true;
}

// =====================================================================
// properOverlap — matches Scala exactly
// =====================================================================
std::string GapFiller::properOverlap(const std::string& left, const std::string& right, int minOverlap) {
    int ll = (int)left.length();
    int rl = (int)right.length();
    for (int overlap = minOverlap; overlap <= ll + rl - 2 * minOverlap; overlap++) {
        int leftOffset = std::max(ll - overlap, 0);
        int rightOffset = std::max(overlap - ll, 0);
        int len = std::min(ll - leftOffset, rl - rightOffset);
        if (substrEq(left, leftOffset, right, rightOffset, len)) {
            return left.substr(0, leftOffset) + right.substr(rightOffset);
        }
    }
    return "";
}

// =====================================================================
// trimPatch — matches Scala exactly
// =====================================================================
std::tuple<int, std::string, std::string>
GapFiller::trimPatch(int startArg, const std::string& patchArg, int stopArg) {
    int start = startArg;
    int stop  = stopArg;
    std::string patch = patchArg;

    while (start < stop && !patch.empty() &&
           (region_.baseAt(start) == patch[0] || region_.originalBaseAt(start) == patch[0])) {
        start++;
        patch = patch.substr(1);
    }
    while (start < stop && !patch.empty() &&
           (region_.baseAt(stop - 1) == patch.back() ||
            region_.originalBaseAt(stop - 1) == patch.back())) {
        stop--;
        patch = patch.substr(0, patch.length() - 1);
    }

    // Build reference string for this range (matching Scala: region.refSubString(start, stop-start))
    int refLen = stop - start;
    std::string ref = region_.refSubString(start, refLen);
    return {start, ref, patch};
}

// =====================================================================
// breakJoins — matches Scala
// =====================================================================
std::vector<std::tuple<int, std::string, std::string>>
GapFiller::breakJoins(int start, const std::vector<std::string>& forwardPaths,
                      const std::vector<std::string>& reversePaths, int stop) {
    std::set<std::tuple<int, std::string, std::string>> solutionSet;

    for (const auto& f : forwardPaths) {
        for (const auto& r : reversePaths) {
            auto s = joinBreak(start, f, r, stop);
            if (s != noSolution) {
                solutionSet.insert(s);
            }
        }
    }

    // Sort by total change length
    std::vector<std::tuple<int, std::string, std::string>> solutions(
        solutionSet.begin(), solutionSet.end());
    std::sort(solutions.begin(), solutions.end(),
              [](const auto& a, const auto& b) {
                  return std::get<2>(a).length() + std::get<1>(a).length()
                       < std::get<2>(b).length() + std::get<1>(b).length();
              });

    // If all solutions have the same delta length, return only the smallest total
    std::set<int> deltas;
    for (const auto& s : solutions) {
        deltas.insert((int)(std::get<2>(s).length()) - (int)(std::get<1>(s).length()));
    }
    if (deltas.size() == 1) {
        return {solutions.front()};
    }
    return solutions;
}

// =====================================================================
// joinBreak — matches Scala
// =====================================================================
std::tuple<int, std::string, std::string>
GapFiller::joinBreak(int startArg, const std::string& forward,
                     const std::string& reverse, int stopArg) {
    int kk = GapFiller::k;
    std::string patch = properOverlap(forward, reverse, kk);

    if (!patch.empty()) {
        auto solution = trimPatch(startArg, patch, stopArg);
        if (std::get<1>(solution) == std::get<2>(solution)) {
            return solution;  // no change, but trimPatch result
        }
        return solution;
    }
    return noSolution;
}

// =====================================================================
// consensusFromLeft — matches Scala
// =====================================================================
std::string GapFiller::consensusFromLeft(const std::vector<std::string>& seqs) const {
    if (seqs.empty()) return "";
    const std::string& s0 = seqs[0];
    if (seqs.size() == 1) return s0;

    int minLength = s0.length();
    for (const auto& s : seqs) minLength = std::min(minLength, (int)s.length());

    for (int i = 0; i < minLength; i++) {
        for (size_t j = 1; j < seqs.size(); j++) {
            if (seqs[j][i] != s0[i]) {
                return s0.substr(0, i);
            }
        }
    }
    return s0.substr(0, minLength);
}

// =====================================================================
// consensusFromRight — matches Scala
// =====================================================================
std::string GapFiller::consensusFromRight(const std::vector<std::string>& seqs) const {
    if (seqs.empty()) return "";
    const std::string& s0 = seqs[0];
    if (seqs.size() == 1) return s0;

    int minLength = s0.length();
    for (const auto& s : seqs) minLength = std::min(minLength, (int)s.length());

    for (int i = 0; i < minLength; i++) {
        for (size_t j = 1; j < seqs.size(); j++) {
            const std::string& sj = seqs[j];
            int si = (int)sj.length() - 1 - i;
            int s0i = (int)s0.length() - 1 - i;
            if (sj[si] != s0[s0i]) {
                return s0.substr(s0i + 1);
            }
        }
    }
    return s0.substr(s0.length() - minLength);
}

// =====================================================================
// partialMatchesReference — matches Scala
// =====================================================================
bool GapFiller::partialMatchesReference(int start, const std::string& fromLeft,
                                        const std::string& fromRight, int stop,
                                        int /*loopLength*/) {
    if ((int)fromLeft.length() > region_.size() || (int)fromRight.length() > region_.size())
        return false;

    if (start + (int)fromLeft.length() > region_.size() ||
        stop - (int)fromRight.length() < 0)
        return false;

    bool leftMatch = true;
    for (int i = 0; i < (int)fromLeft.length() && i + start < region_.size(); i++) {
        if (fromLeft[i] != region_.baseAt(start + i)) { leftMatch = false; break; }
    }
    bool rightMatch = true;
    int rightStart = stop - (int)fromRight.length();
    for (int i = 0; i < (int)fromRight.length() && rightStart + i < region_.size(); i++) {
        if (fromRight[i] != region_.baseAt(rightStart + i)) { rightMatch = false; break; }
    }
    return leftMatch && rightMatch;
}

// =====================================================================
// assembleIntoBreak — matches Scala
// =====================================================================
std::tuple<int, std::vector<std::string>, std::vector<std::string>, int, std::string>
GapFiller::assembleIntoBreak(const Region& brk, const std::vector<BamRead>& reads) {
    Assembler assembler;
    assembler.addReads(reads);
    assembler.buildGraph();

    if (Pilon::fixNovel && !Pilon::novelContigs.empty()) {
        assembler.addGraphSeqs(Pilon::novelContigs);
    }

    int startOffset = breakRadius();
    int start = std::max(region_.start, brk.start - startOffset);
    int stop  = std::min(region_.start + region_.size(), brk.stop + startOffset);

    std::string left  = region_.subString(start, brk.start - start);
    std::string right = region_.subString(brk.stop, stop - brk.stop);

    auto [forward, reverse, loop] = assembler.multiBridge(left, right);
    return {start, forward, reverse, stop, loop};
}

// =====================================================================
// assembleAcrossBreak — matches Scala
// =====================================================================
std::tuple<int, std::string, std::string>
GapFiller::assembleAcrossBreak(const Region& brk, bool isGap) {
    auto reads = recruitReads(brk);
    auto [start, pathsFromLeft, pathsFromRight, stop, loop] = assembleIntoBreak(brk, reads);
    tandemRepeat_ = loop;

    auto solutions = breakJoins(start, pathsFromLeft, pathsFromRight, stop);

    auto solution = noSolution;
    if (solutions.size() == 1 || (Pilon::multiClosure && solutions.size() > 1)) {
        solution = solutions.front();
    }

    bool solutionOK = (solution != noSolution) && (loop.empty() || !Pilon::trSafe);
    if (solutionOK && isGap) {
        int closedLength = (int)std::get<2>(solution).length();
        int closedDiff = std::abs(closedLength - brk.size());
        if (closedDiff > Pilon::gapMargin) solutionOK = false;
    }

    if (solutionOK) {
        return solution;
    }

    if (isGap || (Pilon::fixBreaks && loop.empty())) {
        auto fromRight = consensusFromRight(pathsFromRight);
        auto fromLeft  = consensusFromLeft(pathsFromLeft);

        int newStart = start + (int)fromLeft.length();
        int newStop  = stop  - (int)fromRight.length();

        if ((newStart >= brk.start + GapFiller::minExtend ||
             newStop  <= brk.stop  - GapFiller::minExtend) &&
            !partialMatchesReference(start, fromLeft, fromRight, stop, (int)loop.length())) {

            int newGapLen = isGap ? std::max(Pilon::minGap, newStop - newStart) : Pilon::minGap;
            std::string newGap(newGapLen, 'N');
            std::string seq = fromLeft + newGap + fromRight;
            auto partialSolution = trimPatch(start, seq, stop);
            if (!std::get<2>(partialSolution).empty())
                return partialSolution;
        }
    }
    return noSolution;
}

// =====================================================================
// fillGap / fixBreak — matches Scala
// =====================================================================
std::tuple<int, std::string, std::string> GapFiller::fillGap(const Region& gap) {
    return assembleAcrossBreak(gap, true);
}

std::tuple<int, std::string, std::string> GapFiller::fixBreak(const Region& brk) {
    return assembleAcrossBreak(brk, false);
}

// =====================================================================
// breakRadius — matches Scala
// =====================================================================
int GapFiller::breakRadius() const {
    int minRadius = 3 * Assembler::K;
    double totalInsert = 0.0;
    int count = 0;
    for (auto* bam : BAMS()) {
        if (bam && bam->bamType() == "frags") {
            totalInsert += bam->insertSizeMean();
            count++;
        }
    }
    int insertMean = 0;
    if (count > 0) insertMean = static_cast<int>(std::round(totalInsert / count));
    return std::max(minRadius, insertMean);
}

// =====================================================================
// Read recruitment — matches Scala
// =====================================================================
std::vector<BamRead> GapFiller::recruitReadsFromBams(const Region& reg,
                                                     const std::vector<BamFile*>& bams) const {
    std::vector<BamRead> reads;
    for (auto* b : bams) {
        if (!b) continue;
        auto flankReads = b->recruitFlankReads(reg);
        reads.insert(reads.end(), flankReads.begin(), flankReads.end());
    }
    return reads;
}

std::vector<BamRead> GapFiller::recruitReadsOfType(const Region& reg, const std::string& type) const {
    std::vector<BamFile*> typedBams;
    for (auto* bam : BAMS()) {
        if (bam && bam->bamType() == type) typedBams.push_back(bam);
    }
    return recruitReadsFromBams(reg, typedBams);
}

std::vector<BamRead> GapFiller::recruitFrags(const Region& reg) const {
    return recruitReadsOfType(reg, "frags");
}

std::vector<BamRead> GapFiller::recruitJumps(const Region& reg) const {
    std::vector<BamRead> reads;
    for (auto* bam : BAMS()) {
        if (bam && bam->bamType() == "jumps") {
            auto badMates = bam->recruitBadMates(reg);
            reads.insert(reads.end(), badMates.begin(), badMates.end());
        }
    }
    if (Pilon::debug)
        std::cout << "# Recruiting jump bad mates: count=" << reads.size() << std::endl;
    return reads;
}

std::vector<BamRead> GapFiller::recruitUnpaired(const Region& reg) const {
    if (Pilon::longread) return {};
    return recruitReadsOfType(reg, "unpaired");
}

std::vector<BamRead> GapFiller::recruitReads(const Region& brk) const {
    auto frags = recruitFrags(brk);
    auto jumps = recruitJumps(brk);
    auto unp = recruitUnpaired(brk);
    frags.insert(frags.end(), jumps.begin(), jumps.end());
    frags.insert(frags.end(), unp.begin(), unp.end());
    return frags;
}

// =====================================================================
// Static wrappers (return Fix tuples for fixIssues pipeline)
// =====================================================================
std::tuple<int, std::string, std::string>
GapFiller::doFixGap(GenomeRegion& region, const Region& gap) {
    GapFiller filler(region, region.bamHandles);
    k = 2 * Assembler::K + 1;
    return filler.fillGap(gap);
}

std::tuple<int, std::string, std::string>
GapFiller::doFixBreak(GenomeRegion& region, const Region& brk) {
    GapFiller filler(region, region.bamHandles);
    k = 2 * Assembler::K + 1;
    return filler.fixBreak(brk);
}

// =====================================================================
// closeCircle — matches Scala GapFiller.closeCircle
// =====================================================================
std::vector<std::tuple<int, std::string, std::string>>
GapFiller::closeCircle(int estimatedLength) {
    const int estimatedLengthSlop = 50;
    
    int trimToLength = (estimatedLength > 0) ? estimatedLength : region_.size() / 2;
    int trimFlanks = std::max((region_.size() - trimToLength) / 2, 0);
    int rightEnd = region_.start + trimFlanks;
    Region rightFlank(region_.name, rightEnd, rightEnd + breakRadius() + 1);
    int leftEnd = region_.start + region_.size() - trimFlanks;
    Region leftFlank(region_.name, leftEnd - breakRadius(), leftEnd);
    
    if (Pilon::verbose) {
        std::cout << "left: " << leftFlank.toString() << std::endl;
        std::cout << "right: " << rightFlank.toString() << std::endl;
    }
    
    k = 2 * Assembler::K + 1;
    
    // Recruit reads
    std::vector<BamRead> reads;
    if (estimatedLength == 0) {
        Region fullRegion(region_.name, region_.start, region_.start + region_.size());
        reads = recruitReads(fullRegion);
    } else {
        for (auto* bam : BAMS()) {
            if (bam && bam->bamType() == "unpaired") {
                auto leftReads = bam->readsInRegion(leftFlank);
                auto rightReads = bam->readsInRegion(rightFlank);
                reads.insert(reads.end(), leftReads.begin(), leftReads.end());
                reads.insert(reads.end(), rightReads.begin(), rightReads.end());
            }
        }
    }
    
    if (Pilon::verbose)
        std::cout << "recruited " << reads.size() << " reads" << std::endl;
    
    std::string left = region_.subString(leftFlank.start, leftFlank.size());
    std::string right = region_.subString(rightFlank.start, rightFlank.size());
    
    Assembler assembler;
    assembler.addReads(reads);
    assembler.buildGraph();
    
    auto [forward, reverse, loop] = assembler.multiBridge(left, right);
    if (Pilon::verbose && !loop.empty())
        std::cout << " loop: " << loop.size() << " " << loop << std::endl;
    
    std::set<std::string> patches;
    for (const auto& f : forward) {
        for (const auto& r : reverse) {
            std::string patch = properOverlap(f, r, k);
            if (!patch.empty() &&
                (estimatedLength == 0 ||
                 std::abs(static_cast<int>(patch.size()) - 2 * breakRadius()) < estimatedLengthSlop)) {
                patches.insert(patch);
            }
        }
    }
    
    if (Pilon::verbose) {
        for (const auto& patch : patches)
            std::cout << "patch " << patch.size() << ": " << patch << std::endl;
    }
    
    std::set<int> lengths;
    for (const auto& p : patches) lengths.insert(static_cast<int>(p.size()));
    
    if (Pilon::verbose)
        std::cout << patches.size() << " patches; lengths ";
    for (int l : lengths) std::cout << l << " ";
    if (Pilon::verbose) std::cout << std::endl;
    
    if (lengths.size() == 1 &&
        (estimatedLength == 0 || loop.empty() ||
         std::abs(static_cast<int>(loop.size()) - estimatedLength) < estimatedLengthSlop)) {
        std::string patch = *patches.begin();
        if (Pilon::verbose) std::cout << "  " << patch << std::endl;
        
        // Scala: rightSolution = (1, region.subString(1, rightEnd + breakRadius), "") — 1-based locus
        //        1 in 1-based = contig start, rightEnd = trimFlanks
        // C++: use region_.start as contig start (0-based), same length
        int rightLen = trimFlanks + breakRadius();
        std::vector<std::tuple<int, std::string, std::string>> solutions;
        
        solutions.push_back({region_.start,
                              region_.subString(region_.start, rightLen),
                              ""});
        
        auto leftSolution = trimPatch(leftFlank.start, patch, region_.start + region_.size());
        solutions.push_back(leftSolution);
        
        return solutions;
    }
    
    tandemRepeat_ = loop;
    return {};
}

std::vector<std::tuple<int, std::string, std::string>>
GapFiller::doCloseCircle(GenomeRegion& region, int estimatedLength) {
    GapFiller filler(region, region.bamHandles);
    return filler.closeCircle(estimatedLength);
}

// =====================================================================
// fixNovel — matches Scala + existing implementation
// =====================================================================
void GapFiller::fixNovel(GenomeFile* genome, const std::vector<BamFile*>& /*bamFiles*/) {
    if (Pilon::verbose)
        std::cout << "  fixNovel: assembling novel contigs from unmapped reads" << std::endl;

    Assembler genomeGraph(1);
    if (genome) {
        for (auto& contig : genome->getContigs()) {
            genomeGraph.addGraphSeq(contig.second);
        }
    }

    std::vector<BamRead> unaligned;
    for (auto* bam : Pilon::bamFiles) {
        if (!bam) continue;
        if (bam->bamType() == "jumps") continue;
        auto ua = bam->getUnalignedReads();
        unaligned.insert(unaligned.end(), ua.begin(), ua.end());
    }

    if (unaligned.empty()) return;

    Assembler novelGraph(Pilon::minDepth);
    novelGraph.addReads(unaligned);
    auto novelContigs = novelGraph.novel(genomeGraph);

    // Store in Pilon::novelContigs (used by main.cpp for FASTA output)
    Pilon::novelContigs = novelContigs;
}

} // namespace pilon
