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
#include "pileup.h"

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
    std::string originalBases;
    std::string bases;  // Fixed sequence (output)
    std::vector<PileUp> pileUps;
    double minDepth;

    // Physical coverage tracking
    long long physCovStart;
    long long insertSizeStart;

    // Fix record: position, old sequence, new sequence
    using Fix = std::tuple<int, std::string, std::string>;
    std::vector<Fix> fixes;

    // Read count tracking
    long long readCount;
    long long baseCount;

    // Change tracking (matching Scala changeMap)
    enum ChangeKind { SNP, INS, DEL, AMB };
    struct Change {
        int index;
        ChangeKind kind;
        PileUp pu;
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

    GenomeRegion(const std::string& name, int start, int stop,
                 const std::string& bases, double minDepth);

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

    // Core fixing logic
    void identifyAndFixIssues();

    // Write methods
    void writeVcf(FILE* writer) const;
    void writeChanges(FILE* writer, const std::string& newName, int& offset) const;
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
