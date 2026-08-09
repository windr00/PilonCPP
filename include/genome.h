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

#ifndef PILON_GENOME_H
#define PILON_GENOME_H

#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <cstdint>
#include <functional>
#include <memory>
#include "pileup.h"
#include "normal.h"

namespace pilon {

// Forward declaration
class BamFile;

// Represents a genomic region (gap, break, etc.)
struct Region {
    std::string name;
    int start;  // 0-based
    int stop;   // 0-based exclusive

    Region() : start(0), stop(0) {}
    Region(const std::string& name, int start, int stop)
        : name(name), start(start), stop(stop) {}

    int size() const { return stop - start; }
    int midpoint() const { return start + size() / 2; }

    std::string toString() const {
        std::ostringstream oss;
        oss << name << ":" << start << "-" << stop;
        return oss.str();
    }
};

// Represents a genomic region with pileup data
struct GenomeRegion {
    std::string name;
    int start;       // 0-based start
    int stop;        // 0-based stop (exclusive)
    std::string contigBases;
    std::string fullContig_;   // Full contig sequence for GC sliding window (matching Scala)
    std::string originalBases;
    std::string bases;  // Fixed sequence (output)
    std::vector<PileUp> pileUps;
    double minDepth;

    // Thread-local BAM handles for GapFiller (set before identifyAndFixIssues)
    std::vector<BamFile*>* bamHandles = nullptr;

    // Physical coverage tracking
    long long physCovStart;
    long long insertSizeStart;

    // Fix record: position, old sequence, new sequence
    using Fix = std::tuple<int, std::string, std::string>;
    std::vector<Fix> fixes;
    // Classified fix lists for VCF writeFixRecord filtering (matching Scala)
    std::vector<Fix> bigFixList_;
    std::vector<Fix> smallFixList_;
    // SNP fix list for VCF writeFixRecord filtering (matching Scala snpFixList)
    std::vector<Fix> snpFixList_;
    // Duplication events for VCF <DUP> records (matching Scala)
    std::vector<Region> duplicationEvents_;
    // Reassembly fixes for BED track (region -> outcome label)
    std::vector<std::pair<Region, std::string>> reassemblyFixes_;

    // Read count tracking
    long long readCount;
    long long baseCount;

    // Change tracking (matching Scala changeMap)
    // Stores pre-computed BaseCall data instead of copying entire PileUp
    enum ChangeKind { SNP, INS, DEL, AMB };
    struct Change {
        int index;
        ChangeKind kind;
        // Pre-computed from BaseCall, avoids copying ~180-byte PileUp
        char base;
        bool homo;
        int score;
        char altBase;
        std::string insertion;
        std::string deletion;
        bool homoIndel;
        bool indel;
    };
    std::vector<Change> changes_;

    // Disposition arrays (matching Scala)
    std::vector<bool> confirmed;
    std::vector<bool> ambiguous;
    std::vector<bool> changed;
    std::vector<bool> deleted;
    std::vector<bool> excluded;

    // Summary stats arrays (matching Scala)
    std::vector<int> coverage_arr;
    std::vector<int> badCoverage_arr;
    std::vector<int> clips_arr;
    std::vector<int> insertSize_arr;
    std::vector<int> physCoverage_arr;
    std::vector<int> fragCoverage_arr;
    std::vector<int8_t> weightedQual_arr;
    std::vector<int8_t> weightedMq_arr;
    std::vector<double> copyNumber_arr;   // copy number estimate
    std::vector<int8_t> gc_arr;           // GC content (0-100)

    // GC sliding window buffer
    std::vector<int8_t> gcBuffer;
    int gcCount;

    // Normal distribution objects (matching Scala lazy vals)
    std::unique_ptr<NormalDistribution> physCoverageDist;
    std::unique_ptr<NormalDistribution> coverageDist;
    std::unique_ptr<NormalDistribution> fragCoverageDist;
    std::unique_ptr<NormalDistribution> badCoverageDist;
    std::unique_ptr<NormalDistribution> insertSizeDist;
    std::unique_ptr<NormalDistribution> weightedMqDist;

    // pctBadOverall (matching Scala lazy val)
    int pctBadOverall_;

    GenomeRegion(const std::string& name, int start, int stop,
                 const std::string& bases, const std::string& fullContig,
                 double minDepth);

    int size() const;
    char baseAt(int pos) const;
    char originalBaseAt(int pos) const;
    char refBase(int pos) const;
    std::string subString(int start, int length) const;
    std::string refSubString(int start, int length) const;
    PileUp& pileUpRegion(int index);
    const PileUp& pileUpRegion(int index) const;
    int locus(int index) const;

    // Region helpers (matching Scala Region class)
    bool inRegion(int locus) const;
    bool beforeRegion(int locus) const;
    bool afterRegion(int locus) const;
    int index(int locus) const;

    // Homo run and nanopore exclude (matching Scala)
    int homoRun(int loc) const;
    bool nanoporeExclude(int idx) const;
    void excludeMotifs();

    // Post-processing
    void computePhysCov();
    void postProcess();
    void computeGc(int window = 100);
    void computeCopyNumber();

    // Track-style output (matching Scala Tracks)
    void writeWiggle(FILE* writer, const std::string& name, const std::string& desc,
                     std::function<int(int)> valueFn, const std::string& extraOpts = "") const;
    void writeBed(FILE* writer, const std::string& name,
                  std::function<bool(int)> selectFn) const;
    void writeTracks() const;

    // Accessors for tracks
    int cov(int i) const { return i >= 0 && i < (int)coverage_arr.size() ? coverage_arr[i] : 0; }
    int badCov(int i) const { return i >= 0 && i < (int)badCoverage_arr.size() ? badCoverage_arr[i] : 0; }
    int physCov(int i) const { return i >= 0 && i < (int)physCoverage_arr.size() ? physCoverage_arr[i] : 0; }
    int clip(int i) const { return i >= 0 && i < (int)clips_arr.size() ? clips_arr[i] : 0; }
    int insSize(int i) const { return i >= 0 && i < (int)insertSize_arr.size() ? insertSize_arr[i] : 0; }
    int fragCov(int i) const { return i >= 0 && i < (int)fragCoverage_arr.size() ? fragCoverage_arr[i] : 0; }
    int8_t wQual(int i) const { return i >= 0 && i < (int)weightedQual_arr.size() ? weightedQual_arr[i] : 0; }
    int8_t wMq(int i) const { return i >= 0 && i < (int)weightedMq_arr.size() ? weightedMq_arr[i] : 0; }
    double copyNum(int i) const { return i >= 0 && i < (int)copyNumber_arr.size() ? copyNumber_arr[i] : 0.0; }
    int8_t gc(int i) const { return i >= 0 && i < (int)gc_arr.size() ? gc_arr[i] : 0; }
    bool isChanged(int i) const { return i >= 0 && i < (int)changed.size() ? changed[i] : false; }
    bool isAmbiguous(int i) const { return i >= 0 && i < (int)ambiguous.size() ? ambiguous[i] : false; }
    bool isConfirmed(int i) const { return i >= 0 && i < (int)confirmed.size() ? confirmed[i] : true; }
    bool isDeleted(int i) const { return i >= 0 && i < (int)deleted.size() ? deleted[i] : false; }

    // Core fixing logic
    void identifyAndFixIssues();
    void fixBreakRegion(const Region& breakRegion, const std::string& patch);

    // Break detection helpers (matching Scala)
    int pctBadOverall() const;
    bool lowCoverage(int i) const;
    bool highClipping(int i) const;
    bool tooBad(int i) const;
    double dipFraction(int i) const;
    double dipCoverage(int i, int radius = 100) const;
    int delta(int i, const std::vector<int>& values, int radius = 100) const;
    int deltaCoverage(int i, int radius = 100) const;
    double deltaFraction(int i) const;
    bool breakp(int i) const;
    bool nearEdge(const Region& r, int radius = 100) const;
    bool nearAny(const std::vector<Region>& regions, int distance = 100) const;
    std::vector<Region> summaryRegions(std::function<bool(int)> positionTest, int slop = 100) const;
    std::vector<Region> gaps() const;

    // BED classification methods (matching Scala)
    std::vector<Region> changeRegions() const;
    std::vector<Region> unConfirmedRegions() const;
    std::vector<Region> lowCoverageRegions() const;
    std::vector<Region> possibleCollapsedRepeats() const;
    std::vector<Region> possibleBreaks() const;
    // Write a list of regions as BED with label and RGB color
    void regionsToBed(FILE* writer, const std::vector<Region>& regions,
                      const std::string& label, const std::string& rgb) const;

    // Write methods
    void writeVcf(class Vcf* vcf);
    void writeChanges(FILE* writer, const std::string& newName, int& offset) const;
    
    // Free per-chunk memory after processing (pileups are not needed downstream)
    void freeMemory() {
        pileUps.clear();
        pileUps.shrink_to_fit();
        confirmed.clear();
        ambiguous.clear();
        changed.clear();
        deleted.clear();
        excluded.clear();
        coverage_arr.clear();
        badCoverage_arr.clear();
        clips_arr.clear();
        insertSize_arr.clear();
        physCoverage_arr.clear();
        fragCoverage_arr.clear();
        weightedQual_arr.clear();
        weightedMq_arr.clear();
        homoRunLen_.clear();
        // Full contig sequence is only needed for the GC sliding window during
        // postProcess; it is dead weight once per-chunk processing finishes.
        fullContig_.clear();
        fullContig_.shrink_to_fit();
    }

private:
    // Lazily-precomputed homo-run lengths (homoRun len at each position of
    // contigBases); computed on first homoRun() call, matching the linear scans.
    mutable std::vector<int> homoRunLen_;
};

// FASTA genome file reader
class GenomeFile {
public:
    GenomeFile(const std::string& path, const std::string& targets = "");

    // Returns list of {name, sequence} pairs
    std::vector<std::pair<std::string, std::string>> getContigs() const;

    // Returns list of {name, length} pairs for VCF header
    std::vector<std::pair<std::string, int>> getContigSizes() const;

    void processRegions(std::vector<BamFile*>& bamFiles);

    // Get processed regions for VCF output
    std::vector<GenomeRegion>& getProcessedRegions() { return processedRegions_; }
    const std::vector<GenomeRegion>& getProcessedRegions() const { return processedRegions_; }

private:
    std::vector<std::pair<std::string, std::string>> contigs_;
    std::string targets_;
    std::vector<GenomeRegion> processedRegions_;

    void parseFasta(const std::string& path);
    std::vector<std::pair<std::string, std::string>> parseTargets();
};

} // namespace pilon

#endif // PILON_GENOME_H
