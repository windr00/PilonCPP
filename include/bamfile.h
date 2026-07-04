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

#ifndef PILON_BAMFILE_H
#define PILON_BAMFILE_H

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <hts.h>
#include <sam.h>
#include "assembler.h"
#include "genome.h"

namespace pilon {

class BamFile {
public:
    static constexpr const char * indexSuffix = ".bai";
    static constexpr int minOrientationPct = 10;
    static constexpr int maxFragInsertSize = 700;

    // Long read type codes
    static constexpr int notLongRead = 0;
    static constexpr int nanoporeLongRead = 1;
    static constexpr int pacbioLongRead = 2;

    // Max insert sizes by type
    static std::unordered_map<std::string, int> getMaxInsertSizes();

    BamFile(const std::string& path, const std::string& bamType,
            const std::string& subType = "none");

    ~BamFile();

    const std::string& path() const { return path_; }
    const std::string& bamType() const { return bamType_; }
    const std::string& subType() const { return subType_; }
    int longReadType() const { return longReadType_; }

    // Open and validate the BAM file
    bool open();
    void close();

    // Get sequence names from BAM header
    std::unordered_set<std::string> getSeqNames() const;

    // Process reads in a region
    double process(GenomeRegion& region, int printInterval = 100000);

    // Scan entire BAM for statistics
    void scan(const std::unordered_set<std::string>& seqsOfInterest);

    // Get reads in a region
    std::vector<BamRead> readsInRegion(const Region& region) const;

    // Recruit reads flanking a region
    std::vector<BamRead> recruitFlankReads(const Region& region) const;

    // Recruit bad/stray mates for jump library (matching Scala)
    std::vector<BamRead> recruitBadMates(const Region& region) const;

    // Get unmapped reads for novel contig assembly
    std::vector<BamRead> getUnalignedReads() const;

    // Insert size statistics
    double insertSizeMean() const;
    double insertSizeSigma() const;
    int maxInsertSize() const;

    // Auto-detect BAM type
    std::string autoBam() const;

    // Read counts
    long long mapped() const { return mapped_; }
    long long unmapped() const { return unmapped_; }
    long long filtered() const { return filtered_; }
    long long proper() const { return proper_; }

private:
    std::string path_;
    std::string bamType_;
    std::string subType_;
    int longReadType_;
    long long baseCount_;

    // htslib handles
    htsFile* htsFile_;
    bam_hdr_t* header_;
    hts_idx_t* index_;

    // Statistics
    long long mapped_;
    long long unmapped_;
    long long filtered_;
    long long secondary_;
    long long proper_;

    // Insert size stats
    struct InsertSizeStats {
        long long count = 0;
        long long sum = 0;
        long long sumSq = 0;

        void add(int size);
        double mean() const;
        double sigma() const;
        void reset();
        int maxInsertSize() const;
        std::string toString() const;
    };

    InsertSizeStats insertStatsFR_;
    InsertSizeStats insertStatsRF_;
    InsertSizeStats insertStatsUnpaired_;

    bool validateRead(const bam1_t* read) const;
    void addInsert(int insertSize, bool rc = false, bool unpaired = false);
    int pctFR() const;
    int pctRF() const;
    int pctUnpaired() const;

    Region flankRegion(const Region& region) const;

    // Stray mate tracking
    struct StrayMateMap {
        std::unordered_map<std::string, BamRead> readMap1; // first of pair
        std::unordered_map<std::string, BamRead> readMap2; // second of pair

        void addRead(const BamRead& read) {
            if (read.firstOfPair) {
                readMap1[read.readName] = read;
            } else {
                readMap2[read.readName] = read;
            }
        }

        // Find mates that are missing from the other map
        std::vector<BamRead> findStrays() const {
            std::vector<BamRead> strays;
            for (const auto& [name, read] : readMap1) {
                if (readMap2.find(name) == readMap2.end()) {
                    strays.push_back(read);
                }
            }
            for (const auto& [name, read] : readMap2) {
                if (readMap1.find(name) == readMap1.end()) {
                    strays.push_back(read);
                }
            }
            return strays;
        }

        // Look up mate for a given read
        BamRead* lookupMate(const BamRead& read) const {
            const auto& otherMap = read.firstOfPair ? readMap2 : readMap1;
            auto it = otherMap.find(read.readName);
            if (it != otherMap.end()) return const_cast<BamRead*>(&it->second);
            return nullptr;
        }

        int nStrays() const {
            int count = 0;
            for (const auto& [name, read] : readMap1) {
                if (readMap2.find(name) != readMap2.end()) count++;
            }
            return count * 2;
        }
    };

    mutable StrayMateMap strayMateMap_;
};

} // namespace pilon

#endif // PILON_BAMFILE_H
