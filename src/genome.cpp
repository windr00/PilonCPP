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
#include "genome.h"
#include "utils.h"
#include "pilon.h"
#include "vcf.h"
#include "bases.h"
#include "gapfiller.h"
#include "bamfile.h"
#include "gapfiller.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cctype>
#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>
#include <queue>
#include <functional>
#include <vector>
#include <algorithm>
#include <cmath>
namespace pilon {
// Thread-safe output mutex
static std::mutex coutMutex;
// =============================================================================
// Helper: process a single chunk (called from worker threads)
// =============================================================================
static void processChunk(const std::string& name,
                         const std::string& seq,
                         int chunkStart, int chunkStop,
                         std::vector<BamFile*>& bamFiles,
                         std::vector<GenomeRegion>& results,
                         std::mutex& resultsMutex,
                         std::atomic<int>& completedChunks,
                         int totalChunks) {
    std::vector<BamFile*> threadBams;
    for (size_t i = 0; i < bamFiles.size(); i++) {
        auto* bam = bamFiles[i];
        if (bam) {
            auto* threadBam = new BamFile(bam->path(), bam->bamType(), bam->subType());
            threadBam->open();
            threadBam->shareScanState(*bam);
            threadBams.push_back(threadBam);
        }
    }
    GenomeRegion region(name, chunkStart, chunkStop,
                        seq.substr(chunkStart, chunkStop - chunkStart),
                        seq,  // full contig for GC sliding window
                        Pilon::minDepth);
    for (auto* bam : threadBams) {
        std::vector<long long> covBefore;
        if (bam->bamType() != "jumps") {
            covBefore.reserve(region.size());
            for (int i = 0; i < region.size(); i++)
                covBefore.push_back(region.pileUpRegion(i).depth());
        }
        bam->process(region);
        if (bam->bamType() != "jumps") {
            for (int i = 0; i < region.size(); i++)
                region.fragCoverage_arr[i] += static_cast<int>(region.pileUpRegion(i).depth() - covBefore[i]);
        }
    }
    region.bamHandles = &threadBams;
    region.postProcess();
    if (Pilon::fixSnps || Pilon::fixIndels || Pilon::fixGaps || Pilon::fixLocal) {
        region.identifyAndFixIssues();
    }
    {
        std::lock_guard<std::mutex> lock(resultsMutex);
        results.push_back(std::move(region));
    }
    int done = ++completedChunks;
    if (Pilon::verbose || done % 10 == 0 || done == totalChunks) {
        std::lock_guard<std::mutex> lock(coutMutex);
        std::cout << "  [" << done << "/" << totalChunks << "] Chunk "
                  << chunkStart << "-" << chunkStop << " done" << std::endl;
    }
    for (auto* bam : threadBams) {
        bam->close();
        delete bam;
    }
}
// =============================================================================
// GenomeRegion constructor
// =============================================================================
GenomeRegion::GenomeRegion(const std::string& name, int start, int stop,
                           const std::string& bases, const std::string& fullContig,
                           double minDepth)
    : name(name), start(start), stop(stop), contigBases(bases),
      fullContig_(fullContig), originalBases(bases), bases(bases), minDepth(minDepth),
      physCovStart(0), insertSizeStart(0), readCount(0), baseCount(0),
      pctBadOverall_(0) {
    int sz = stop - start;
    pileUps.resize(sz);
    confirmed.resize(sz, false);
    ambiguous.resize(sz, false);
    changed.resize(sz, false);
    deleted.resize(sz, false);
    excluded.resize(sz, false);
    coverage_arr.resize(sz, 0);
    badCoverage_arr.resize(sz, 0);
    clips_arr.resize(sz, 0);
    insertSize_arr.resize(sz, 0);
    physCoverage_arr.resize(sz, 0);
    fragCoverage_arr.resize(sz, 0);
    weightedQual_arr.resize(sz, 0);
    weightedMq_arr.resize(sz, 0);
}
// =============================================================================
// Basic accessors
// =============================================================================
int GenomeRegion::size() const { return stop - start; }
char GenomeRegion::baseAt(int pos) const {
    int local = pos - start;
    // Scala: defaults to mutable bases (post-SNP-correction) for subString/trimPatch
    if (local >= 0 && local < static_cast<int>(bases.size()))
        return bases[local];
    return 'N';
}
char GenomeRegion::originalBaseAt(int pos) const {
    int local = pos - start;
    if (local >= 0 && local < static_cast<int>(originalBases.size()))
        return originalBases[local];
    return 'N';
}
char GenomeRegion::refBase(int pos) const {
    // Scala: returns from contigBases (immutable original reference), not mutable bases
    int local = pos - start;
    if (local >= 0 && local < static_cast<int>(contigBases.size()))
        return contigBases[local];
    return 'N';
}
std::string GenomeRegion::subString(int locus, int length) const {
    int local = locus - start;
    // Scala: uses mutable bases (updated by SNP fixes before gap/break assembly)
    if (local < 0 || local >= static_cast<int>(bases.size())) return "";
    int len = std::min(length, static_cast<int>(bases.size()) - local);
    return bases.substr(local, len);
}
std::string GenomeRegion::refSubString(int locus, int length) const {
    int local = locus - start;
    // Scala: uses immutable contigBases (original reference)
    if (local < 0 || local >= static_cast<int>(contigBases.size())) return "";
    int len = std::min(length, static_cast<int>(contigBases.size()) - local);
    return contigBases.substr(local, len);
}
PileUp& GenomeRegion::pileUpRegion(int index) {
    return pileUps[index];
}
const PileUp& GenomeRegion::pileUpRegion(int index) const {
    return pileUps[index];
}
int GenomeRegion::locus(int index) const {
    return start + index;
}
bool GenomeRegion::inRegion(int locus) const {
    return locus >= start && locus < stop;
}
bool GenomeRegion::beforeRegion(int locus) const {
    return locus < start;
}
bool GenomeRegion::afterRegion(int locus) const {
    return locus >= stop;
}
int GenomeRegion::index(int locus) const {
    return locus - start;
}
// =============================================================================
// Homo run and nanopore exclude (matching Scala)
// =============================================================================
int GenomeRegion::homoRun(int loc) const {
    if (loc < 0 || loc >= static_cast<int>(contigBases.size())) return 0;
    // One-pass precompute of run length at every position (reverse scan).
    if (homoRunLen_.empty()) {
        int n = static_cast<int>(contigBases.size());
        homoRunLen_.assign(n, 1);
        for (int i = n - 2; i >= 0; i--) {
            homoRunLen_[i] = (contigBases[i] == contigBases[i + 1]) ? homoRunLen_[i + 1] + 1 : 1;
        }
    }
    return homoRunLen_[loc];
}
bool GenomeRegion::nanoporeExclude(int idx) const {
    return (idx - 2 >= 0 && idx + 2 < static_cast<int>(contigBases.size()) &&
            contigBases[idx - 2] == 'C' &&
            contigBases[idx - 1] == 'C' &&
            contigBases[idx + 1] == 'G' &&
            contigBases[idx + 2] == 'G');
}
void GenomeRegion::excludeMotifs() {
    bool nano = Pilon::nanopore;
    for (int i = 0; i < static_cast<int>(contigBases.size()); i++) {
        excluded[i] = (homoRun(i) >= 4) || (nano && nanoporeExclude(i));
    }
}
// =============================================================================
// Physical coverage computation (matching Scala PileUpRegion.computePhysCov)
// =============================================================================
void GenomeRegion::computePhysCov() {
    if (pileUps.empty()) return;
    
    pileUps[0].physCov += physCovStart;
    pileUps[0].insertSize += insertSizeStart;
    
    for (int i = 1; i < static_cast<int>(pileUps.size()); i++) {
        pileUps[i].physCov += pileUps[i - 1].physCov;
        pileUps[i].insertSize += pileUps[i - 1].insertSize;
    }
    
    for (int i = 0; i < static_cast<int>(pileUps.size()); i++) {
        if (pileUps[i].physCov > 0) {
            pileUps[i].insertSize /= pileUps[i].physCov;
        }
    }
}
// =============================================================================
// GC content computation (matching Scala sliding window)
// =============================================================================
void GenomeRegion::computeGc(int window) {
    gc_arr.resize(size());
    gcBuffer.resize(window, 0);

    // Initialize circular buffer to ~50% GC
    gcCount = 0;
    for (int i = 0; i < window; i++) {
        if ((i & 1) == 0) {
            gcBuffer[i] = 1;
            gcCount++;
        } else {
            gcBuffer[i] = 0;
        }
    }

    int halfWindow = (window + 1) / 2;
    // Use full contig for sliding window (matching Scala: for (locus <- 0 until contig.length))
    for (int locus = 0; locus < (int)fullContig_.size(); locus++) {
        int center = (locus >= halfWindow) ? locus - halfWindow : locus;
        int bufIndex = locus % window;
        char base = fullContig_[locus];
        int8_t gcBase;
        if (base == 'G' || base == 'C') gcBase = 1;
        else if (base == 'A' || base == 'T') gcBase = 0;
        else gcBase = gcBuffer[bufIndex]; // NOP for Ns, IUPAC, etc

        gcCount += gcBase - gcBuffer[bufIndex];
        gcBuffer[bufIndex] = gcBase;

        if (inRegion(center))
            gc_arr[index(center)] = static_cast<int8_t>(gcCount);

        if (inRegion(locus) && locus >= static_cast<int>(fullContig_.size()) - halfWindow)
            gc_arr[index(locus)] = static_cast<int8_t>(gcCount);
    }
}

// =============================================================================
// Copy number estimation (matching Scala)
// =============================================================================
// Smooth function: sliding window average (matching Scala smooth())
static std::vector<int> smooth(const std::vector<int>& input, int window) {
    int inputSize = static_cast<int>(input.size());
    std::vector<int> result(inputSize, 0);
    int half = window / 2;
    long long accum = 0;
    for (int i = 0; i < inputSize; i++) {
        accum += input[i];
        if (i > window) {
            accum -= input[i - window];
            result[i - half] = (accum + half) / window;
        }
    }
    if (inputSize > window) {
        for (int i = 0; i < window - half; i++) result[i] = result[window - half];
        for (int i = inputSize - half; i < inputSize; i++) result[i] = result[inputSize - half - 1];
    } else {
        for (int i = 0; i < inputSize; i++) result[i] = accum / inputSize;
    }
    return result;
}

void GenomeRegion::computeCopyNumber() {
    copyNumber_arr.resize(size(), 1.0);
    if (size() == 0) return;
    auto smoothCov = smooth(fragCoverage_arr, 200);
    double baseCov = fragCoverageDist ? fragCoverageDist->mean : 1.0;
    for (int i = 0; i < size(); i++) {
        double cn = baseCov > 0 ? smoothCov[i] / baseCov : 0.0;
        copyNumber_arr[i] = std::round(cn);
    }
}

// =============================================================================
// BED classification methods (matching Scala Tracks + GenomeRegion)
// =============================================================================

void GenomeRegion::regionsToBed(FILE* writer, const std::vector<Region>& regions,
                                const std::string& label, const std::string& rgb) const {
    if (!writer) return;
    for (const auto& region : regions) {
        fprintf(writer, "%s\t%d\t%d\t%s\t0\t+\t%d\t%d\t%s\n",
                region.name.c_str(), region.start - 1, region.stop,
                label.c_str(), region.start - 1, region.stop, rgb.c_str());
    }
}

std::vector<Region> GenomeRegion::changeRegions() const {
    return summaryRegions([this](int i) -> bool { return isChanged(i) || isDeleted(i); }, 1);
}

std::vector<Region> GenomeRegion::unConfirmedRegions() const {
    return summaryRegions([this](int i) -> bool { return !isConfirmed(i); });
}

std::vector<Region> GenomeRegion::lowCoverageRegions() const {
    return summaryRegions([this](int i) -> bool { return lowCoverage(i); });
}

std::vector<Region> GenomeRegion::possibleCollapsedRepeats() const {
    auto highCopy = summaryRegions([this](int i) -> bool { return copyNum(i) > 1; });
    std::vector<Region> result;
    for (const auto& r : highCopy) {
        int startIdx = index(r.start);
        int stopIdx = index(r.stop);
        if (startIdx >= 0 && stopIdx >= 0 && deltaFraction(startIdx) >= 0.5 && deltaFraction(stopIdx) >= 0.5) {
            result.push_back(r);
        }
    }
    return result;
}

std::vector<Region> GenomeRegion::possibleBreaks() const {
    auto breaks = summaryRegions([this](int i) -> bool { return breakp(i); }, 200);
    auto gps = gaps();
    std::vector<Region> result;
    for (const auto& brk : breaks) {
        // Check if break is near any gap (matching Scala: filter { !_.nearAny(gaps, 300) })
        bool nearGap = false;
        for (const auto& g : gps) {
            int distance = std::abs(brk.start - g.stop);
            distance = std::min(distance, std::abs(brk.stop - g.start));
            if (distance <= 300) { nearGap = true; break; }
        }
        if (!nearGap) result.push_back(brk);
    }
    return result;
}

// =============================================================================
// Track output (matching Scala Tracks class)
// =============================================================================
void GenomeRegion::writeWiggle(FILE* writer, const std::string& name,
                               const std::string& desc,
                               std::function<int(int)> valueFn,
                               const std::string& extraOpts) const {
    if (!writer) return;
    fprintf(writer, "track type=wiggle_0 name=\"%s\" description=\"%s\" %s\n",
            name.c_str(), desc.c_str(), extraOpts.c_str());
    fprintf(writer, "fixedStep chrom=%s start=%d step=1\n", name.c_str(), start + 1);
    for (int i = 0; i < size(); i++) {
        fprintf(writer, "%d\n", valueFn(i));
    }
}

void GenomeRegion::writeBed(FILE* writer, const std::string& name,
                            std::function<bool(int)> selectFn) const {
    if (!writer) return;
    int bedStart = -1;
    for (int i = 0; i < size(); i++) {
        if (selectFn(i)) {
            if (bedStart < 0) bedStart = i;
        } else {
            if (bedStart >= 0) {
                fprintf(writer, "%s\t%d\t%d\t%s\n",
                        name.c_str(), start + bedStart, start + i, name.c_str());
                bedStart = -1;
            }
        }
    }
    if (bedStart >= 0) {
        fprintf(writer, "%s\t%d\t%d\t%s\n",
                name.c_str(), start + bedStart, start + size(), name.c_str());
    }
}

void GenomeRegion::writeTracks() const {
    // BED track (matching Scala: 8 color-classified region tracks)
    std::string bedFile = Pilon::outputFile(".bed");
    FILE* bed = fopen(bedFile.c_str(), "a");
    if (bed) {
        // Write BED header on first contig (check file position)
        if (ftell(bed) == 0) {
            fprintf(bed, "track description=\"Issues found by Pilon\" name=\"Pilon\"\n");
        }
        regionsToBed(bed, unConfirmedRegions(), "?", "0,128,128");
        regionsToBed(bed, changeRegions(), "X", "255,0,0");
        regionsToBed(bed, lowCoverageRegions(), "LowCov", "128,0,128");
        regionsToBed(bed, possibleCollapsedRepeats(), "Copy#", "128,128,0");
        regionsToBed(bed, duplicationEvents_, "Duplication", "255,128,0");
        regionsToBed(bed, gaps(), "Gap", "0,0,0");
        regionsToBed(bed, possibleBreaks(), "Break", "255,0,255");
        for (const auto& rf : reassemblyFixes_) {
            std::vector<Region> single = {rf.first};
            regionsToBed(bed, single, rf.second, "0,128,0");
        }
        fclose(bed);
    }

    // Changes.wig
    std::string chFn = Pilon::outputFile("Changes.wig");
    FILE* chF = fopen(chFn.c_str(), "a");
    if (chF) {
        writeWiggle(chF, name, "Changes",
                    [this](int i) { return isChanged(i) ? 1 : 0; });
        fclose(chF);
    }

    // Unconfirmed.wig
    std::string ucFn = Pilon::outputFile("Unconfirmed.wig");
    FILE* ucF = fopen(ucFn.c_str(), "a");
    if (ucF) {
        writeWiggle(ucF, name, "Unconfirmed",
                    [this](int i) { return isConfirmed(i) ? 0 : 1; });
        fclose(ucF);
    }

    // CopyNumber.wig
    std::string cnFn = Pilon::outputFile("CopyNumber.wig");
    FILE* cnF = fopen(cnFn.c_str(), "a");
    if (cnF) {
        writeWiggle(cnF, name, "Copy Number",
                    [this](int i) { return static_cast<int>(copyNum(i) - 1); });
        fclose(cnF);
    }

    // Coverage.wig
    std::string covFn = Pilon::outputFile("Coverage.wig");
    FILE* covF = fopen(covFn.c_str(), "a");
    if (covF) {
        writeWiggle(covF, name, "Coverage", [this](int i) { return cov(i); });
        fclose(covF);
    }

    // BadCoverage.wig
    std::string bcFn = Pilon::outputFile("BadCoverage.wig");
    FILE* bcF = fopen(bcFn.c_str(), "a");
    if (bcF) {
        writeWiggle(bcF, name, "Bad Coverage", [this](int i) { return badCov(i); });
        fclose(bcF);
    }

    // PctBad.wig
    std::string pbFn = Pilon::outputFile("PctBad.wig");
    FILE* pbF = fopen(pbFn.c_str(), "a");
    if (pbF) {
        writeWiggle(pbF, name, "Pct Bad",
                    [this](int i) {
                        int good = cov(i);
                        int bad = badCov(i);
                        return (good + bad > 0) ? bad * 100 / (good + bad) : 0;
                    });
        fclose(pbF);
    }

    // DeltaCoverage.wig
    std::string dcFn = Pilon::outputFile("DeltaCoverage.wig");
    FILE* dcF = fopen(dcFn.c_str(), "a");
    if (dcF) {
        writeWiggle(dcF, name, "Delta Coverage",
                    [this](int i) { return deltaCoverage(i, 100); });
        fclose(dcF);
    }

    // DipCoverage.wig
    std::string dipFn = Pilon::outputFile("DipCoverage.wig");
    FILE* dipF = fopen(dipFn.c_str(), "a");
    if (dipF) {
        writeWiggle(dipF, name, "Dip Coverage",
                    [this](int i) { return static_cast<int>(dipCoverage(i, 100)); });
        fclose(dipF);
    }

    // PhysicalCoverage.wig
    std::string physFn = Pilon::outputFile("PhysicalCoverage.wig");
    FILE* physF = fopen(physFn.c_str(), "a");
    if (physF) {
        writeWiggle(physF, name, "Physical Coverage", [this](int i) { return physCov(i); });
        fclose(physF);
    }

    // ClippedAlignments.wig
    std::string clipFn = Pilon::outputFile("ClippedAlignments.wig");
    FILE* clipF = fopen(clipFn.c_str(), "a");
    if (clipF) {
        writeWiggle(clipF, name, "Clipped Alignments", [this](int i) { return clip(i); });
        fclose(clipF);
    }

    // WeightedQual.wig
    std::string qFn = Pilon::outputFile("WeightedQual.wig");
    FILE* qF = fopen(qFn.c_str(), "a");
    if (qF) {
        writeWiggle(qF, name, "Weighted Qual", [this](int i) { return wQual(i); });
        fclose(qF);
    }

    // WeightedMq.wig
    std::string mqFn = Pilon::outputFile("WeightedMq.wig");
    FILE* mqF = fopen(mqFn.c_str(), "a");
    if (mqF) {
        writeWiggle(mqF, name, "Weighted MQ", [this](int i) { return wMq(i); });
        fclose(mqF);
    }

    // GC.wig
    std::string gcFn = Pilon::outputFile("GC.wig");
    FILE* gcF = fopen(gcFn.c_str(), "a");
    if (gcF) {
        writeWiggle(gcF, name, "GC",
                    [this](int i) { return gc(i); },
                    "graphType=heatmap midRange=35:65 midColor=0,255,0");
        fclose(gcF);
    }
}

    // =============================================================================
    // fixBreakRegion: apply a local assembly patch to the contig bases
    // =============================================================================
    void GenomeRegion::fixBreakRegion(const Region& breakRegion, const std::string& patch) {
        int startPos = start + breakRegion.start;

        if (Pilon::verbose) {
            std::cout << "    fixing break " << breakRegion.start << "-" << breakRegion.stop
                      << " patch_len=" << patch.length() << std::endl;
        }

        if (startPos > (int)contigBases.size()) return;

        std::string newBases;
        newBases.append(contigBases.substr(0, startPos));
        newBases.append(patch);
        int afterBreak = startPos + (breakRegion.stop - breakRegion.start);
        if (afterBreak < (int)contigBases.size())
            newBases.append(contigBases.substr(afterBreak));

        contigBases = std::move(newBases);

        // Also update bases for output (FASTA writes from bases)
        // Rebuild bases from the new contigBases
        bases = contigBases;

        if (Pilon::verbose)
            std::cout << "    break fixed, new length: " << contigBases.size() << std::endl;
    }
void GenomeRegion::postProcess() {
    computePhysCov();
    
    // Compute meanCoverage (matching Scala: roundDiv(baseCount, size))
    int meanCoverage = 0;
    if (!pileUps.empty()) {
        long long totalBaseCount = 0;
        for (const auto& pu : pileUps) {
            totalBaseCount += pu.baseCount.sum();
        }
        long long d = static_cast<long long>(pileUps.size());
        meanCoverage = static_cast<int>((totalBaseCount + d / 2) / d);
    }
    
    // Compute minDepth (matching Scala postProcess logic)
    minDepth = Pilon::minDepth >= 1 
        ? static_cast<int>(Pilon::minDepth)
        : std::max(static_cast<int>(std::round(Pilon::minDepth * meanCoverage)), Pilon::minMinDepth);
    
    if (readCount == 0) {
        if (Pilon::verbose) {
            std::cerr << "WARNING: No reads processed in region " << name 
                      << ":" << start << "-" << stop << std::endl;
        }
        computeGc();
        computeCopyNumber();
        // Initialize normal distributions with empty data to avoid null deref
        std::vector<int> empty;
        coverageDist.reset(new NormalDistribution(empty, 2));
        fragCoverageDist.reset(new NormalDistribution(empty, 2));
        pctBadOverall_ = 0;
        return;
    }
    
    // excluded motifs for long reads (matching Scala: longReadOnly = Pilon.longread && fragBams.isEmpty)
    {
        bool hasFrags = false;
        for (auto* bam : Pilon::bamFiles) {
            if (bam && bam->bamType() == "frags") { hasFrags = true; break; }
        }
        if (Pilon::longread && !hasFrags) {
            excludeMotifs();
        }
    }
    
    // Pass 1: pull out values from pileups & call base changes
    bool fixamb = Pilon::iupac || Pilon::fixAmb;
    for (int i = 0; i < size(); i++) {
        const PileUp& pu = pileUps[i];
        long long n = pu.depth();
        auto bc = pu.baseCall();
        char b = bc.base;
        bool homo = bc.homo;
        char r = refBase(locus(i));
        
        // Store summary stats (matching Scala)
        coverage_arr[i] = n;
        badCoverage_arr[i] = pu.badPair;
        physCoverage_arr[i] = pu.physCov;
        insertSize_arr[i] = pu.insertSize;
        weightedQual_arr[i] = static_cast<int8_t>(pu.weightedQual());
        weightedMq_arr[i] = static_cast<int8_t>(pu.weightedMq());
        clips_arr[i] = pu.clips;
        
        if (n >= minDepth && r != 'N' && !deleted[i] && bc.called()) {
            if (homo && b == r && bc.highConfidence() && !bc.indel) {
                confirmed[i] = true;
            } else if (bc.isInsertion() && bc.homoIndel) {
                changed[i] = true;
                // Compact Change: pre-compute BaseCall data, no PileUp copy
                changes_.push_back({i, INS, b, homo, bc.score,
                                    static_cast<char>(0), bc.insertion, bc.deletion,
                                    bc.homoIndel, bc.indel});
            } else if (bc.isDeletion() && bc.homoIndel) {
                changed[i] = true;
                // Compact Change: pre-compute BaseCall data, no PileUp copy
                changes_.push_back({i, DEL, b, homo, bc.score,
                                    static_cast<char>(0), bc.insertion, bc.deletion,
                                    bc.homoIndel, bc.indel});
                // Scala: mark deleted positions for downstream positions
                for (int j = 1; j < static_cast<int>(bc.deletion.length()); j++) {
                    int delIdx = i + j;
                    if (delIdx < static_cast<int>(deleted.size())) {
                        deleted[delIdx] = true;
                        pileUps[delIdx].deletions += pu.deletions;
                    }
                }
            } else if (b != r && bc.score > 0) {
                if (homo) {
                    changed[i] = true;
                    // SNP Change
                    changes_.push_back({i, SNP, b, homo, bc.score,
                                        bc.altBase, std::string(), std::string(),
                                        bc.homoIndel, bc.indel});
                } else if (fixamb || bc.altBase != r) {
                    ambiguous[i] = true;
                    // AMB Change
                    changes_.push_back({i, AMB, b, homo, bc.score,
                                        bc.altBase, std::string(), std::string(),
                                        bc.homoIndel, bc.indel});
                }
            } else {
            }
        }
    }
    
    // Pass 2: computed values
    computeGc();

    // Compute normal distributions (matching Scala lazy vals)
    if (!coverage_arr.empty()) {
        coverageDist.reset(new NormalDistribution(coverage_arr, 2));
        badCoverageDist.reset(new NormalDistribution(badCoverage_arr, 2));
        physCoverageDist.reset(new NormalDistribution(physCoverage_arr, 2));
        insertSizeDist.reset(new NormalDistribution(insertSize_arr, 2));
        weightedMqDist.reset(new NormalDistribution(
            std::vector<int>(weightedMq_arr.begin(), weightedMq_arr.end()), 2));

        // fragCoverage: use non-jump coverage (deep copy then compute)
        std::vector<int> fragVec;
        fragVec.reserve(fragCoverage_arr.size());
        for (int v : fragCoverage_arr) fragVec.push_back(v);
        fragCoverageDist.reset(new NormalDistribution(fragVec, 2));
    }

    // Copy number estimation (requires fragCoverageDist to be initialized)
    computeCopyNumber();

    // Compute pctBadOverall (matching Scala lazy val)
    long long totalBad = 0;
    long long totalGood = 0;
    for (int i = 0; i < size(); i++) {
        totalBad += badCoverage_arr[i];
        totalGood += coverage_arr[i];
    }
    pctBadOverall_ = Utils::pct(totalBad, totalGood + totalBad);
}
// =============================================================================
// fixFixList: Sort and remove overlapping fixes (matching Scala)
// =============================================================================
static std::vector<GenomeRegion::Fix> fixFixList(const std::vector<GenomeRegion::Fix>& inList) {
    std::vector<GenomeRegion::Fix> fixes = inList;
    std::sort(fixes.begin(), fixes.end(),
              [](const GenomeRegion::Fix& a, const GenomeRegion::Fix& b) {
                  return std::get<0>(a) < std::get<0>(b);
              });
    
    std::vector<GenomeRegion::Fix> outList;
    
    while (!fixes.empty()) {
        if (fixes.size() >= 2) {
            const auto& fix1 = fixes[0];
            const auto& fix2 = fixes[1];
            
            int region1End = std::get<0>(fix1) + std::max(static_cast<int>(std::get<1>(fix1).length()) - 1, 0);
            int region2End = std::get<0>(fix2) + std::max(static_cast<int>(std::get<1>(fix2).length()) - 1, 0);
            
            if (std::get<0>(fix1) <= region2End && std::get<0>(fix2) <= region1End) {
                int fix1len = std::get<1>(fix1).length() + std::get<2>(fix1).length();
                int fix2len = std::get<1>(fix2).length() + std::get<2>(fix2).length();
                
                if (fix1len >= fix2len) {
                    // Scala: fixes = fix1 :: tail  (keep fix1, drop fix2, compare fix1 with tail.head next)
                    fixes.erase(fixes.begin() + 1);
                } else {
                    // Scala: fixes = fix2 :: tail  (keep fix2, drop fix1, compare fix2 with tail.head next)
                    fixes.erase(fixes.begin());
                }
            } else {
                outList.push_back(fix1);
                fixes.erase(fixes.begin());
            }
        } else {
            outList.push_back(fixes[0]);
            fixes.clear();
        }
    }
    
    return outList;
}

// =============================================================================
// Break detection helpers (matching Scala)
// =============================================================================

int GenomeRegion::pctBadOverall() const { return pctBadOverall_; }

int GenomeRegion::delta(int i, const std::vector<int>& values, int radius) const {
    int left = values[std::max(0, i - radius)];
    int right = values[std::min(size() - 1, i + radius)];
    return std::abs(left - right);
}

int GenomeRegion::deltaCoverage(int i, int radius) const {
    return delta(i, fragCoverage_arr, radius);
}

double GenomeRegion::deltaFraction(int i) const {
    double mean = fragCoverageDist ? fragCoverageDist->mean : 1.0;
    if (mean <= 0) return 0;
    return deltaCoverage(i) / mean;
}

bool GenomeRegion::lowCoverage(int i) const {
    return coverage_arr[i] < Pilon::minMinDepth && refBase(locus(i)) != 'N';
}

bool GenomeRegion::highClipping(int i) const {
    return coverage_arr[i] >= Pilon::minMinDepth &&
           Utils::pct(static_cast<long long>(clips_arr[i]),
                      static_cast<long long>(coverage_arr[i])) >= 33;
}

bool GenomeRegion::tooBad(int i) const {
    int good = coverage_arr[i];
    int bad = badCoverage_arr[i];
    int p = Utils::pct(static_cast<long long>(bad),
                       static_cast<long long>(good + bad));
    return p > pctBadOverall_ + 20;
}

double GenomeRegion::dipCoverage(int i, int radius) const {
    int left = fragCoverage_arr[std::max(0, i - radius)];
    int right = fragCoverage_arr[std::min(size() - 1, i + radius)];
    int center = fragCoverage_arr[i];
    return (left - center) + (right - center);
}

double GenomeRegion::dipFraction(int i) const {
    double mean = fragCoverageDist ? fragCoverageDist->mean : 1.0;
    if (mean <= 0) return 0;
    return dipCoverage(i) / mean;
}

bool GenomeRegion::breakp(int i) const {
    return lowCoverage(i) || highClipping(i) || tooBad(i) || (dipFraction(i) >= 1.5);
}

bool GenomeRegion::nearEdge(const Region& r, int radius) const {
    return r.start - start < radius || stop - r.stop < radius;
}

bool GenomeRegion::nearAny(const std::vector<Region>& regions, int distance) const {
    // nearAny is called per-break from possibleBreaks; use region bounds to compare
    // However, since this is called from possibleBreaks() which iterates individual breaks,
    // we need to check each gap against the chunk. To fix: inline the check in possibleBreaks.
    for (const auto& other : regions) {
        if (other.name == name) {
            // Use chunk bounds as a fallback check (gap proximity to chunk)
            // Individual break checks are handled inline in possibleBreaks()
            if (std::abs(other.stop - start) <= distance ||
                std::abs(other.start - stop) <= distance) return true;
        }
    }
    return false;
}

std::vector<Region> GenomeRegion::summaryRegions(std::function<bool(int)> positionTest, int slop) const {
    std::vector<Region> regions;
    int first = -1;
    int last = -1;

    for (int i = 0; i < size(); i++) {
        if (positionTest(i)) {
            last = i;
            if (first < 0) first = i;
        } else {
        if (last >= 0 && i > last + slop) {
            regions.push_back(Region(name, locus(first), locus(last) + 1));
            first = -1;
            last = -1;
        }
    }
}
if (last >= 0) regions.push_back(Region(name, locus(first), locus(last) + 1));

    // Reverse and filter near-edge (matching Scala's regions.reverse.filter(!nearEdge(_)))
    std::reverse(regions.begin(), regions.end());
    std::vector<Region> result;
    for (const auto& r : regions) {
        if (!nearEdge(r)) result.push_back(r);
    }

    return result;
}

std::vector<Region> GenomeRegion::gaps() const {
    auto gapTest = [this](int i) -> bool {
        return refBase(locus(i)) == 'N';
    };
    auto rawGaps = summaryRegions(gapTest);
    std::vector<Region> result;
    for (const auto& g : rawGaps) {
        if (g.size() >= 10) result.push_back(g);
    }
    return result;
}

// =============================================================================
// fixIssues: Apply fixes to bases (matching Scala fixIssues)
// =============================================================================
void GenomeRegion::identifyAndFixIssues() {
    std::vector<Fix> snpFixList;
    std::vector<Fix> smallFixList;
    std::vector<Fix> bigFixList;
    
    if (Pilon::verbose) {
        std::cout << "  Changes found: " << changes_.size() << std::endl;
    }
    
    for (const auto& change : changes_) {
        int i = change.index;
        ChangeKind kind = change.kind;
        
        int loc = locus(i);
        char rBase = refBase(loc);
        char cBase = change.base;
        
        if (!excluded[i]) {
            switch (kind) {
                case SNP:
                    if (Pilon::fixSnps) {
                        snpFixList.push_back({loc, std::string(1, rBase), std::string(1, cBase)});
                    }
                    break;
                case AMB:
                    if (Pilon::fixSnps && !Pilon::longread) {
                        if (Pilon::iupac) {
                            // Scala: Bases.toIUPAC(cBase, bc.altBase)
                            char iupacB = Bases::toIUPAC(cBase, change.altBase);
                            smallFixList.push_back({loc, std::string(1, rBase), std::string(1, iupacB)});
                        } else {
                            snpFixList.push_back({loc, std::string(1, rBase), std::string(1, cBase)});
                        }
                    }
                    break;
                case INS:
                    if (Pilon::fixIndels) {
                        smallFixList.push_back({loc, "", change.insertion});
                    }
                    break;
                case DEL:
                    if (Pilon::fixIndels) {
                        smallFixList.push_back({loc, change.deletion, ""});
                    }
                    break;
            }
        }
    }
    
    // Report stats (matching Scala)
    int ins = 0, dels = 0;
    int insBases = 0, delBases = 0;
    for (const auto& fix : smallFixList) {
        if (std::get<1>(fix).empty() && !std::get<2>(fix).empty()) {
            ins++;
            insBases += std::get<2>(fix).length();
        } else if (!std::get<1>(fix).empty() && std::get<2>(fix).empty()) {
            dels++;
            delBases += std::get<1>(fix).length();
        }
    }
    
    // Apply SNP fixes first (matching Scala: fixIssues(snpFixList))
    // Do this BEFORE gap/break assembly so they work on corrected bases
    if (!snpFixList.empty()) {
        auto snpFixes = fixFixList(snpFixList);
        for (auto it = snpFixes.rbegin(); it != snpFixes.rend(); ++it) {
            int locus_global = std::get<0>(*it);
            const std::string& was = std::get<1>(*it);
            const std::string& patch = std::get<2>(*it);
            int startIdx = locus_global - start;
            for (size_t i = 0; i < was.length(); i++) {
                int iNew = startIdx + i;
                if (iNew >= 0 && iNew < static_cast<int>(bases.length())) {
                    bases[iNew] = patch[i];
                }
            }
        }
    }

    // Try to fill gaps (matching Scala identifyAndFixIssues)
    if (Pilon::fixGaps) {
        auto gapRegions = gaps();
        if (!gapRegions.empty()) {
            if (Pilon::compact) {
                printf("  # Filling %zu gaps: ...", gapRegions.size());
                fflush(stdout);
            } else {
                printf("  # Filling %zu gaps:\n", gapRegions.size());
            }
            int gapIdx = 0;
            for (const auto& gap : gapRegions) {
                gapIdx++;
                bool showDetail = (Pilon::verbose && !Pilon::compact) || (!Pilon::compact && (gapIdx % 5 == 0 || gapIdx == 1));
                if (Pilon::compact && gapIdx % 5 == 0) {
                    printf("\r  # Filling gaps: %d/%zu  ", gapIdx, gapRegions.size());
                    fflush(stdout);
                }
                if (showDetail) {
                    printf("    gap %d/%zu: %s", gapIdx, gapRegions.size(), gap.toString().c_str());
                }
                auto fix = GapFiller::doFixGap(*this, gap);
                int fixStart = std::get<0>(fix);
                const std::string& ref = std::get<1>(fix);
                const std::string& patch = std::get<2>(fix);
                if (fixStart > 0) {
                    bigFixList.push_back(fix);
                    if (showDetail)
                        printf(" -> %zubp patch", patch.length());
                } else {
                    if (showDetail)
                        printf(" -> no solution");
                }
                if (showDetail)
                    printf("\n");
                // Populate reassemblyFixes_ (matching Scala logFix)
                {
                    int nRef = 0, nPatch = 0;
                    for (char c : ref) if (c == 'N') nRef++;
                    for (char c : patch) if (c == 'N') nPatch++;
                    std::string msg;
                    if (fixStart == 0) msg = "NoSolution";
                    else if (nRef == gap.size() && nPatch == 0) msg = "ClosedGap";
                    else if (gap.size() > 0) msg = "PartialFill";
                    else msg = "Unknown!";
                    reassemblyFixes_.push_back({gap, msg});
                }
            }
            if (Pilon::compact) {
                printf("\r  # Filling gaps: %zu/%zu done\n", gapRegions.size(), gapRegions.size());
            }
        }
    }

    // Try to fix local breaks (matching Scala identifyAndFixIssues)
    if (Pilon::fixLocal) {
        auto breakRegions = summaryRegions([this](int i) -> bool { return breakp(i); }, 200);
        // Filter breaks near gaps (matching Scala: filter { !_.nearAny(gaps, 300) })
        auto gapRegions = gaps();
        std::vector<Region> filteredBreaks;
        for (const auto& brk : breakRegions) {
            // Check nearAny against gapRegions
            bool nearGap = false;
            for (const auto& g : gapRegions) {
                if (std::abs(g.start - brk.stop) <= 300 || std::abs(g.stop - brk.start) <= 300 ||
                    (brk.start <= g.stop && g.start <= brk.stop)) {
                    nearGap = true;
                    break;
                }
            }
            if (!nearGap) filteredBreaks.push_back(brk);
        }
        if (!filteredBreaks.empty()) {
            if (Pilon::compact) {
                printf("  # Fixing %zu breaks: ...", filteredBreaks.size());
                fflush(stdout);
            } else {
                printf("  # Fixing %zu breaks:\n", filteredBreaks.size());
            }
            int brkIdx = 0;
            for (const auto& brk : filteredBreaks) {
                brkIdx++;
                bool showDetail = (Pilon::verbose && !Pilon::compact) || (!Pilon::compact && (brkIdx % 5 == 0 || brkIdx == 1));
                if (Pilon::compact && brkIdx % 5 == 0) {
                    printf("\r  # Fixing breaks: %d/%zu  ", brkIdx, filteredBreaks.size());
                    fflush(stdout);
                }
                if (showDetail)
                    printf("    break %d/%zu: %s", brkIdx, filteredBreaks.size(), brk.toString().c_str());
                auto fix = GapFiller::doFixBreak(*this, brk);
                int fixStart = std::get<0>(fix);
                const std::string& ref = std::get<1>(fix);
                const std::string& patch = std::get<2>(fix);
                if (fixStart > 0 && std::max(static_cast<int>(ref.length()), static_cast<int>(patch.length())) > 10) {
                    bigFixList.push_back(fix);
                    if (showDetail)
                        printf(" -> %zubp patch", patch.length());
                } else {
                    if (showDetail)
                        printf(" -> no solution");
                }
                if (showDetail)
                    printf("\n");
                // Populate reassemblyFixes_ (matching Scala logFix)
                {
                    int nRef = 0, nPatch = 0;
                    for (char c : ref) if (c == 'N') nRef++;
                    for (char c : patch) if (c == 'N') nPatch++;
                    std::string msg;
                    if (fixStart == 0) msg = "NoSolution";
                    else if (ref.empty() && patch.empty()) msg = "NoChange";
                    else if (nPatch == 0) msg = "BreakFix";
                    else if (nPatch > 0 && nRef == 0) msg = "OpenedGap";
                    else msg = "Unknown!";
                    reassemblyFixes_.push_back({brk, msg});
                }
            }
            if (Pilon::compact) {
                printf("\r  # Fixing breaks: %zu/%zu done\n", filteredBreaks.size(), filteredBreaks.size());
            }
        }
    }

    // fixCircles: close circular contigs (matching Scala closeCircle)
    // Only apply to first chunk (contig origin); Scala region starts at contig base 1
    if (Pilon::fixCircles && start == 0) {
        auto circleFixes = GapFiller::doCloseCircle(*this, 0);
        for (const auto& fix : circleFixes) {
            bigFixList.push_back(fix);
        }
    }

    // Report large collapsed regions (matching Scala duplicationEvents)
    auto dups = summaryRegions([this](int i) -> bool { return copyNumber_arr[i] > 1; }, 2000);
    // Store filtered duplications (matching Scala: regions filter {_.size > 10000})
    duplicationEvents_.clear();
    for (const auto& d : dups) {
        if (d.size() > 10000) duplicationEvents_.push_back(d);
    }
    if (!dups.empty()) {
        for (const auto& d : dups) {
            if (d.size() > 10000 && Pilon::verbose)
                std::cout << "Large collapsed region: " << d.toString() << " size " << d.size() << std::endl;
        }
    }

    fixes.insert(fixes.end(), snpFixList.begin(), snpFixList.end());
    fixes.insert(fixes.end(), smallFixList.begin(), smallFixList.end());
    fixes.insert(fixes.end(), bigFixList.begin(), bigFixList.end());

    // Store classified fix lists for VCF writeFixRecord and BED filtering (matching Scala)
    snpFixList_ = snpFixList;
    bigFixList_ = bigFixList;
    smallFixList_ = smallFixList;

    // Deduplicate overlapping fixes (matching Scala fixFixList for VCF/changes output)
    fixes = fixFixList(fixes);

    // Apply all fixes (small + big) to bases (matching Scala: fixIssues(smallFixList ++ bigFixList))
    {
        std::vector<Fix> allFixes;
        allFixes.insert(allFixes.end(), smallFixList.begin(), smallFixList.end());
        allFixes.insert(allFixes.end(), bigFixList.begin(), bigFixList.end());
        
        auto sortedFixes = fixFixList(allFixes);
        for (auto it = sortedFixes.rbegin(); it != sortedFixes.rend(); ++it) {
            int locus_global = std::get<0>(*it);
            const std::string& was = std::get<1>(*it);
            const std::string& patch = std::get<2>(*it);
            int startIdx = locus_global - start;
            
            if (was.length() == patch.length()) {
                for (size_t i = 0; i < was.length(); i++) {
                    int iNew = startIdx + i;
                    if (iNew >= 0 && iNew < static_cast<int>(bases.length())) {
                        bases[iNew] = patch[i];
                    }
                }
            } else {
                if (startIdx >= 0 && startIdx <= static_cast<int>(bases.length())) {
                    std::string before = bases.substr(0, startIdx);
                    std::string after = (startIdx + static_cast<int>(was.length()) <= static_cast<int>(bases.length()))
                        ? bases.substr(startIdx + was.length())
                        : "";
                    bases = before + patch + after;
                }
            }
        }
    }
}
// =============================================================================
// Write VCF records (matching Scala GenomeRegion.writeVcf)
// =============================================================================
void GenomeRegion::writeVcf(Vcf* vcf) {
    if (!vcf) return;
    
    // Merge and deduplicate fix lists (matching Scala: fixFixList(snpFixList ++ smallFixList ++ bigFixList))
    // Note: snpFixList entries are already in `fixes`; we need all classified fixes here
    std::vector<Fix> allClassified;
    allClassified.insert(allClassified.end(), snpFixList_.begin(), snpFixList_.end());
    allClassified.insert(allClassified.end(), smallFixList_.begin(), smallFixList_.end());
    allClassified.insert(allClassified.end(), bigFixList_.begin(), bigFixList_.end());
    auto fixesSorted = fixFixList(allClassified);
    size_t fixIdx = 0;
    
    // Duplication events sorted by start position
    // duplicationEvents_ are already stored as Region objects with global coordinates
    auto dupes = duplicationEvents_;
    size_t dupIdx = 0;
    
    for (int i = 0; i < size(); i++) {
        int loc = locus(i);
        
        // Write DUP record if a duplication event starts at this position
        while (dupIdx < dupes.size() && dupes[dupIdx].start < loc) dupIdx++;
        if (dupIdx < dupes.size() && dupes[dupIdx].start == loc) {
            vcf->writeDup(*this, dupes[dupIdx]);
            dupIdx++;
        }
        
        // Write fix record for big fixes (reassembly) not also in small list
        if (fixIdx < fixesSorted.size() && std::get<0>(fixesSorted[fixIdx]) == loc) {
            const auto& fix = fixesSorted[fixIdx];
            bool inBig = false;
            for (const auto& bf : bigFixList_) {
                if (bf == fix) { inBig = true; break; }
            }
            bool inSmall = false;
            for (const auto& sf : smallFixList_) {
                if (sf == fix) { inSmall = true; break; }
            }
            if (inBig && !inSmall) {
                vcf->writeFixRecord(*this, fix);
                for (int j = 0; j < static_cast<int>(std::get<1>(fix).length()); j++) {
                    if (i + j < size()) deleted[i + j] = true;
                }
            }
            fixIdx++;
        }
        
        // Write per-position record
        vcf->writeRecord(*this, i, isDeleted(i));
    }
}

// =============================================================================
// Write changes
// =============================================================================
void GenomeRegion::writeChanges(FILE* writer, const std::string& newName, int& offset) const {
    int delta = 0;
    for (const auto& fix : fixes) {
        int loc = std::get<0>(fix) + offset + 1; // 1-based
        int newLoc = loc + delta;
        const std::string& refSeq = std::get<1>(fix);
        const std::string& altSeq = std::get<2>(fix);
        std::string fromDisplay = refSeq.empty() ? "." : refSeq;
        std::string toDisplay = altSeq.empty() ? "." : altSeq;
        fprintf(writer, "%s\t%d\t%s\t%s\n",
                newName.c_str(), newLoc, fromDisplay.c_str(), toDisplay.c_str());
        delta += static_cast<int>(altSeq.length()) - static_cast<int>(refSeq.length());
    }
    offset += static_cast<int>(bases.size()) - size();
}
// =============================================================================
// GenomeFile implementation
// =============================================================================
GenomeFile::GenomeFile(const std::string& path, const std::string& targets)
    : targets_(targets) {
    parseFasta(path);
}
void GenomeFile::parseFasta(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open genome file: " << path << std::endl;
        exit(1);
    }
    std::string line;
    std::string currentName;
    std::string currentSeq;
    while (std::getline(file, line)) {
        line = Utils::trim(line);
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!currentName.empty()) {
                contigs_.push_back({currentName, currentSeq});
                currentSeq.clear();
            }
            std::istringstream iss(line.substr(1));
            iss >> currentName;
        } else {
            for (char c : line) {
                currentSeq += std::toupper(c);
            }
        }
    }
    if (!currentName.empty()) {
        contigs_.push_back({currentName, currentSeq});
    }
}
std::vector<std::pair<std::string, std::string>> GenomeFile::getContigs() const {
    return contigs_;
}
std::vector<std::pair<std::string, int>> GenomeFile::getContigSizes() const {
    std::vector<std::pair<std::string, int>> sizes;
    for (const auto& contig : contigs_) {
        sizes.push_back({contig.first, static_cast<int>(contig.second.size())});
    }
    return sizes;
}
std::vector<std::pair<std::string, std::string>> GenomeFile::parseTargets() {
    if (targets_.empty()) return contigs_;
    std::string content = Utils::readFile(targets_);
    std::istringstream stream(content);
    std::vector<std::pair<std::string, std::string>> result;
    std::string line;
    while (std::getline(stream, line)) {
        line = Utils::trim(line);
        if (line.empty() || line[0] == '#') continue;
        auto parts = Utils::split(line, '\t');
        if (parts.size() >= 1) {
            std::string name = Utils::trim(parts[0]);
            for (const auto& contig : contigs_) {
                if (contig.first == name) {
                    result.push_back(contig);
                    break;
                }
            }
        }
    }
    return result;
}
void GenomeFile::processRegions(std::vector<BamFile*>& bamFiles) {
    auto contigs = targets_.empty() ? contigs_ : parseTargets();
    int numThreads = Pilon::threads;
    if (numThreads <= 1) {
        for (const auto& contig : contigs) {
            const std::string& name = contig.first;
            const std::string& seq = contig.second;
            int length = static_cast<int>(seq.size());
            std::cout << "Processing " << name << " (" << length << " bp)" << std::endl;
            int numChunks = (length + Pilon::chunkSize - 1) / Pilon::chunkSize;
            int chunkIdx = 0;
            for (int chunkStart = 0; chunkStart < length; chunkStart += Pilon::chunkSize) {
                chunkIdx++;
                int chunkStop = std::min(chunkStart + Pilon::chunkSize, length);
                GenomeRegion region(name, chunkStart, chunkStop,
                                   seq.substr(chunkStart, chunkStop - chunkStart),
                                   seq,  // full contig for GC sliding window
                                   Pilon::minDepth);
                for (auto* bam : bamFiles) {
                    if (bam) {
                        std::vector<long long> covBefore;
                        if (bam->bamType() != "jumps") {
                            covBefore.reserve(region.size());
                            for (int i = 0; i < region.size(); i++)
                                covBefore.push_back(region.pileUpRegion(i).depth());
                        }
                        bam->process(region);
                        if (bam->bamType() != "jumps") {
                            for (int i = 0; i < region.size(); i++)
                                region.fragCoverage_arr[i] += static_cast<int>(region.pileUpRegion(i).depth() - covBefore[i]);
                        }
                    }
                }
                region.bamHandles = &bamFiles;
                region.postProcess();
                if (Pilon::fixSnps || Pilon::fixIndels || Pilon::fixGaps || Pilon::fixLocal) {
                    region.identifyAndFixIssues();
                }
                processedRegions_.push_back(std::move(region));
                // Free per-chunk memory unless VCF output is needed
                if (!Pilon::vcf) {
                    processedRegions_.back().freeMemory();
                }
                if (Pilon::verbose || chunkIdx % 5 == 0 || chunkIdx == numChunks) {
                    std::cout << "  [" << chunkIdx << "/" << numChunks << "] Chunk "
                              << chunkStart << "-" << chunkStop << " done" << std::endl;
                }
            }
        }
        return;
    }
    // Multi-threaded mode
    std::cout << "Using " << numThreads << " threads for processing" << std::endl;
    struct ChunkTask { std::string name; std::string seq; int chunkStart; int chunkStop; };
    std::vector<ChunkTask> tasks;
    for (const auto& contig : contigs) {
        const std::string& name = contig.first;
        const std::string& seq = contig.second;
        int length = static_cast<int>(seq.size());
        std::cout << "Processing " << name << " (" << length << " bp)" << std::endl;
        for (int chunkStart = 0; chunkStart < length; chunkStart += Pilon::chunkSize) {
            int chunkStop = std::min(chunkStart + Pilon::chunkSize, length);
            tasks.push_back({name, seq, chunkStart, chunkStop});
        }
    }
    int totalChunks = static_cast<int>(tasks.size());
    std::cout << "Total chunks: " << totalChunks << std::endl;
    std::mutex resultsMutex;
    std::atomic<int> completedChunks(0);
    processedRegions_.reserve(totalChunks);
    std::vector<std::thread> workers;
    workers.reserve(numThreads);
    std::mutex taskMutex;
    std::queue<size_t> taskQueue;
    for (size_t i = 0; i < tasks.size(); i++) taskQueue.push(i);
    auto workerFunc = [&]() {
        while (true) {
            size_t taskIdx;
            {
                std::lock_guard<std::mutex> lock(taskMutex);
                if (taskQueue.empty()) return;
                taskIdx = taskQueue.front();
                taskQueue.pop();
            }
            const auto& task = tasks[taskIdx];
            processChunk(task.name, task.seq, task.chunkStart, task.chunkStop,
                        bamFiles, processedRegions_, resultsMutex,
                        completedChunks, totalChunks);
        }
    };
    for (int i = 0; i < numThreads; i++) workers.emplace_back(workerFunc);
    for (auto& t : workers) t.join();

    // Natural sort comparator: "scaffold2" < "scaffold10"
    auto naturalLess = [](const std::string& a, const std::string& b) -> bool {
        size_t ai = 0, bi = 0;
        while (ai < a.size() && bi < b.size()) {
            if (std::isdigit(a[ai]) && std::isdigit(b[bi])) {
                // Numeric segment: compare by value
                size_t ae = ai, be = bi;
                while (ae < a.size() && std::isdigit(a[ae])) ae++;
                while (be < b.size() && std::isdigit(b[be])) be++;
                unsigned long long na = std::stoull(a.substr(ai, ae - ai));
                unsigned long long nb = std::stoull(b.substr(bi, be - bi));
                if (na != nb) return na < nb;
                ai = ae; bi = be;
            } else {
                // String segment: lexicographic
                if (a[ai] != b[bi]) return a[ai] < b[bi];
                ai++; bi++;
            }
        }
        return a.size() < b.size();
    };

    std::sort(processedRegions_.begin(), processedRegions_.end(),
              [&](const GenomeRegion& a, const GenomeRegion& b) {
                  if (a.name != b.name) return naturalLess(a.name, b.name);
                  return a.start < b.start;
              });
    
    // Free per-chunk memory unless VCF output is needed
    if (!Pilon::vcf) {
        for (auto& reg : processedRegions_) {
            reg.freeMemory();
        }
    }
    
    std::cout << "Multi-threaded processing complete: " << totalChunks << " chunks" << std::endl;
}
} // namespace pilon
