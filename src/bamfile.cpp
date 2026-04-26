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

#include "bamfile.h"
#include "pilon.h"
#include "utils.h"
#include <iostream>
#include <cmath>
#include <cstring>
#include <algorithm>

namespace pilon {

std::unordered_map<std::string, int> BamFile::getMaxInsertSizes() {
    return {
        {"frags", 500},
        {"jumps", 10000},
        {"unpaired", 10000},
        {"bam", 10000}
    };
}

BamFile::BamFile(const std::string& path, const std::string& bamType,
                 const std::string& subType)
    : path_(path), bamType_(bamType), subType_(subType),
      longReadType_(notLongRead), baseCount_(0),
      htsFile_(nullptr), header_(nullptr), index_(nullptr),
      mapped_(0), unmapped_(0), filtered_(0), secondary_(0), proper_(0) {
    
    if (bamType_ == "unpaired") {
        if (subType_ == "nanopore") longReadType_ = nanoporeLongRead;
        else if (subType_ == "pacbio") longReadType_ = pacbioLongRead;
    }
}

BamFile::~BamFile() {
    close();
}

bool BamFile::open() {
    htsFile_ = hts_open(path_.c_str(), "rb");
    if (!htsFile_) {
        std::cerr << "Error: Cannot open BAM file: " << path_ << std::endl;
        return false;
    }

    header_ = sam_hdr_read(htsFile_);
    if (!header_) {
        std::cerr << "Error: Cannot read BAM header: " << path_ << std::endl;
        hts_close(htsFile_);
        htsFile_ = nullptr;
        return false;
    }

    // Try to open index (.bai or .csi)
    std::string indexPath = path_ + ".bai";
    if (!Utils::fileExists(indexPath)) {
        indexPath = path_ + ".csi";
    }
    
    if (Utils::fileExists(indexPath)) {
        index_ = hts_idx_load2(path_.c_str(), indexPath.c_str());
    }

    return true;
}

void BamFile::close() {
    if (index_) {
        hts_idx_destroy(index_);
        index_ = nullptr;
    }
    if (header_) {
        sam_hdr_destroy(header_);
        header_ = nullptr;
    }
    if (htsFile_) {
        hts_close(htsFile_);
        htsFile_ = nullptr;
    }
}

std::unordered_set<std::string> BamFile::getSeqNames() const {
    std::unordered_set<std::string> names;
    if (!header_) return names;

    // Use sam_hdr_nref for newer htslib, fall back to sam_hdr_nseq
    int nseq = sam_hdr_nref(header_);
    for (int i = 0; i < nseq; i++) {
        const char* name = sam_hdr_tid2name(header_, i);
        if (name) {
            names.insert(name);
        }
    }
    return names;
}

// Helper: check if mate is unmapped using mtid
static bool isMateUnmapped(const bam1_t* bam) {
    // In older htslib, mate_flag doesn't exist; use mtid == -1 as proxy
    return bam->core.mtid < 0;
}

// Helper: convert a bam1_t to BamRead
static BamRead bamToRead(const bam1_t* bam, const bam_hdr_t* hdr) {
    BamRead read;
    read.readName = std::string(bam_get_qname(bam));
    
    // Convert encoded bases (4 bits per base) to ASCII
    std::string decoded;
    decoded.reserve(bam->core.l_qseq);
    const uint8_t* seq = bam_get_seq(bam);
    for (int i = 0; i < bam->core.l_qseq; i++) {
        int base = seq[i] & 0xF;
        if (base == 1) decoded += 'A';
        else if (base == 2) decoded += 'C';
        else if (base == 4) decoded += 'G';
        else if (base == 8) decoded += 'T';
        else decoded += 'N';
    }
    read.bases = decoded;

    const uint8_t* qual = bam_get_qual(bam);
    read.quals.assign(qual, qual + bam->core.l_qseq);
    read.mappingQuality = bam->core.qual;
    read.alignmentStart = bam->core.pos;
    read.alignmentEnd = bam_endpos(bam);
    read.negativeStrand = (bam->core.flag & BAM_FREVERSE) != 0;
    read.paired = (bam->core.flag & BAM_FPAIRED) != 0;
    read.properPair = (bam->core.flag & BAM_FPROPER_PAIR) != 0;
    read.firstOfPair = (bam->core.flag & BAM_FREAD1) != 0;
    read.unmapped = (bam->core.flag & BAM_FUNMAP) != 0;
    read.mateUnmapped = isMateUnmapped(bam);
    read.duplicate = (bam->core.flag & BAM_FDUP) != 0;
    read.secondary = (bam->core.flag & BAM_FSECONDARY) != 0;
    read.failsVendorQC = (bam->core.flag & BAM_FQCFAIL) != 0;
    read.inferredInsertSize = bam->core.isize;

    if (bam->core.tid >= 0) {
        read.referenceName = sam_hdr_tid2name(hdr, bam->core.tid);
    }
    if (bam->core.mtid >= 0) {
        read.mateReferenceName = sam_hdr_tid2name(hdr, bam->core.mtid);
    }
    read.mateAlignmentStart = bam->core.mpos;

    // Parse CIGAR string
    uint32_t ncigar = bam->core.n_cigar;
    const uint32_t* cigar = bam_get_cigar(bam);
    std::string cigarStr;
    for (uint32_t i = 0; i < ncigar; i++) {
        cigarStr += std::to_string(bam_cigar_oplen(cigar[i]));
        cigarStr += static_cast<char>(bam_cigar_opchr(bam_cigar_op(cigar[i])));
    }
    read.cigar = cigarStr;

    return read;
}

// Helper: query BAM for reads in a region using htslib iterator
static std::vector<BamRead> queryRegion(const htsFile* fp, const bam_hdr_t* hdr,
                                         hts_idx_t* idx, const std::string& refName,
                                         int start, int stop) {
    std::vector<BamRead> reads;
    if (!fp || !hdr) return reads;

    // Cast away const for older htslib API
    int tid = sam_hdr_name2tid(const_cast<bam_hdr_t*>(hdr), refName.c_str());
    if (tid < 0) return reads;

    hts_itr_t* iter = nullptr;
    if (idx) {
        iter = sam_itr_queryi(idx, tid, start, stop);
    } else {
        iter = sam_itr_queryi(nullptr, tid, start, stop);
    }
    if (!iter) return reads;

    bam1_t* bam = bam_init1();
    int ret;
    while ((ret = sam_itr_next(const_cast<htsFile*>(fp), iter, bam)) > 0) {
        reads.push_back(bamToRead(bam, hdr));
    }
    if (ret < -1) {
        std::cerr << "Error reading BAM in region " << refName << ":" 
                  << start << "-" << stop << std::endl;
    }

    hts_itr_destroy(iter);
    bam_destroy1(bam);
    return reads;
}

double BamFile::process(const GenomeRegion& region, int printInterval) {
    if (!htsFile_) {
        if (!open()) return 0;
    }

    int tid = sam_hdr_name2tid(header_, region.name.c_str());
    if (tid < 0) return 0;

    hts_itr_t* iter = nullptr;
    if (index_) {
        iter = sam_itr_queryi(index_, tid, region.start, region.stop);
    } else {
        iter = sam_itr_queryi(nullptr, tid, region.start, region.stop);
    }
    if (!iter) return 0;

    bam1_t* bam = bam_init1();
    int lastLoc = 0;

    int ret;
    while ((ret = sam_itr_next(htsFile_, iter, bam)) > 0) {
        int loc = bam->core.pos;
        if (Pilon::verbose && printInterval > 0 && loc > lastLoc + printInterval) {
            lastLoc = printInterval * (loc / printInterval);
            std::cout << "..." << lastLoc << std::flush;
        }

        if (validateRead(bam)) {
            int insertSize = std::abs(bam->core.isize);
            bool rc = (bam->core.flag & BAM_FREVERSE) != 0;
            addInsert(insertSize, rc);
        }
    }

    if (ret < -1) {
        std::cerr << "Error reading BAM file" << std::endl;
    }

    hts_itr_destroy(iter);
    bam_destroy1(bam);

    return 0;
}

void BamFile::scan(const std::unordered_set<std::string>& seqsOfInterest) {
    if (!htsFile_) {
        if (!open()) return;
    }

    bam1_t* bam = bam_init1();
    int ret;

    while ((ret = sam_read1(htsFile_, header_, bam)) > 0) {
        if (!validateRead(bam)) {
            filtered_++;
        } else if (bam->core.flag & BAM_FUNMAP) {
            unmapped_++;
        } else {
            mapped_++;
            if (bam->core.flag & BAM_FPAIRED) {
                if (bam->core.flag & BAM_FPROPER_PAIR) {
                    proper_++;
                    addInsert(std::abs(bam->core.isize),
                             (bam->core.flag & BAM_FREVERSE) != 0);
                } else if (Pilon::strays && !isMateUnmapped(bam)) {
                    const char* refName = sam_hdr_tid2name(header_, bam->core.tid);
                    if (refName && seqsOfInterest.count(refName)) {
                        BamRead read = bamToRead(bam, header_);
                        strayMateMap_.addRead(read);
                    }
                }
            } else {
                addInsert(bam->core.l_qseq,
                         (bam->core.flag & BAM_FREVERSE) != 0, true);
            }
        }
    }

    bam_destroy1(bam);

    long long totalReads = mapped_ + unmapped_ + filtered_;
    std::cout << path_ << ": " << totalReads << " reads, "
              << filtered_ << " filtered, "
              << mapped_ << " mapped, "
              << proper_ << " proper";

    if (Pilon::strays) {
        std::cout << ", " << strayMateMap_.nStrays() << " stray";
    }

    int frPct = pctFR();
    int rfPct = pctRF();
    int unPct = pctUnpaired();

    if (frPct >= minOrientationPct) {
        std::cout << ", FR " << frPct << "% " << insertStatsFR_.toString();
    }
    if (rfPct >= minOrientationPct) {
        std::cout << ", RF " << rfPct << "% " << insertStatsRF_.toString();
    }
    if (unPct >= minOrientationPct) {
        std::cout << ", Unpaired " << unPct << "% " << insertStatsUnpaired_.toString();
    }

    std::cout << ", max " << maxInsertSize() << std::endl;

    if (bamType_ == "bam") {
        bamType_ = autoBam();
        std::cout << "  Auto-detected type: " << bamType_ << std::endl;
    }
}

std::vector<BamRead> BamFile::readsInRegion(const Region& region) const {
    return queryRegion(htsFile_, header_, index_, region.name, region.start, region.stop);
}

std::vector<BamRead> BamFile::recruitFlankReads(const Region& region) const {
    Region flanks = flankRegion(region);
    std::vector<BamRead> reads = readsInRegion(flanks);
    
    // Handle stray mates if Pilon::strays is set
    if (Pilon::strays) {
        // Build a temporary mate map from the flank reads
        std::unordered_map<std::string, BamRead> readMap1;
        std::unordered_map<std::string, BamRead> readMap2;
        for (const auto& r : reads) {
            if (r.firstOfPair) {
                readMap1[r.readName] = r;
            } else {
                readMap2[r.readName] = r;
            }
        }
        
        // Find strays: reads in one map but not the other
        for (const auto& [name, read] : readMap1) {
            if (readMap2.find(name) == readMap2.end()) {
                // Look up mate from the global stray mate map
                BamRead* mate = strayMateMap_.lookupMate(read);
                if (mate) reads.push_back(*mate);
            }
        }
        for (const auto& [name, read] : readMap2) {
            if (readMap1.find(name) == readMap1.end()) {
                BamRead* mate = strayMateMap_.lookupMate(read);
                if (mate) reads.push_back(*mate);
            }
        }
    }
    
    return reads;
}

double BamFile::insertSizeMean() const {
    return pctFR() >= 50 ? insertStatsFR_.mean() : insertStatsRF_.mean();
}

double BamFile::insertSizeSigma() const {
    return pctFR() >= 50 ? insertStatsFR_.sigma() : insertStatsRF_.sigma();
}

int BamFile::maxInsertSize() const {
    int m1 = insertStatsFR_.maxInsertSize();
    int m2 = insertStatsRF_.maxInsertSize();
    int m3 = insertStatsUnpaired_.maxInsertSize();
    int maxFromStats = std::max(m1, std::max(m2, m3));
    
    if (maxFromStats > 0) return maxFromStats;
    
    auto maxSizes = getMaxInsertSizes();
    auto it = maxSizes.find(bamType_);
    if (it != maxSizes.end()) return it->second;
    return 10000;
}

std::string BamFile::autoBam() const {
    int fr = pctFR();
    int rf = pctRF();
    int un = pctUnpaired();

    if (un >= fr && un >= rf) return "unpaired";

    double insertSize = rf > fr ? insertStatsRF_.mean() : insertStatsFR_.mean();
    if (insertSize >= maxFragInsertSize) return "jumps";
    return "frags";
}

bool BamFile::validateRead(const bam1_t* bam) const {
    if (!Pilon::nonPf && (bam->core.flag & BAM_FQCFAIL)) return false;
    if (!Pilon::duplicates && (bam->core.flag & BAM_FDUP)) return false;
    if (bam->core.flag & BAM_FSECONDARY) return false;
    return true;
}

void BamFile::addInsert(int insertSize, bool rc, bool unpaired) {
    const int huge = 5 * 10000;
    
    if (insertSize <= 0 || insertSize >= huge) return;

    bool fr = (insertSize > 0) ^ rc;

    if (unpaired) {
        insertStatsUnpaired_.add(insertSize);
    } else if (fr) {
        insertStatsFR_.add(insertSize);
    } else {
        insertStatsRF_.add(insertSize);
    }
}

int BamFile::pctFR() const {
    long long total = insertStatsFR_.count + insertStatsRF_.count + insertStatsUnpaired_.count;
    if (total == 0) return 0;
    return (insertStatsFR_.count * 100) / total;
}

int BamFile::pctRF() const {
    long long total = insertStatsFR_.count + insertStatsRF_.count + insertStatsUnpaired_.count;
    if (total == 0) return 0;
    return (insertStatsRF_.count * 100) / total;
}

int BamFile::pctUnpaired() const {
    long long total = insertStatsFR_.count + insertStatsRF_.count + insertStatsUnpaired_.count;
    if (total == 0) return 0;
    return (insertStatsUnpaired_.count * 100) / total;
}

Region BamFile::flankRegion(const Region& region) const {
    int flank = maxInsertSize();
    return Region(region.name,
                  std::max(region.start - flank, 0),
                  region.stop + flank);
}

// InsertSizeStats implementation
void BamFile::InsertSizeStats::add(int size) {
    sum += size;
    sumSq += static_cast<long long>(size) * size;
    count++;
}

double BamFile::InsertSizeStats::mean() const {
    if (count == 0) return 0.0;
    return static_cast<double>(sum) / count;
}

double BamFile::InsertSizeStats::sigma() const {
    if (count == 0 || sum == 0) return 0.0;
    double m = mean();
    double variance = (static_cast<double>(sumSq) / count) - (m * m);
    return std::sqrt(std::abs(variance));
}

void BamFile::InsertSizeStats::reset() {
    count = 0;
    sum = 0;
    sumSq = 0;
}

int BamFile::InsertSizeStats::maxInsertSize() const {
    if (count < 1000) return 0;
    return static_cast<int>(mean() + 3.0 * sigma());
}

std::string BamFile::InsertSizeStats::toString() const {
    char buf[64];
    snprintf(buf, sizeof(buf), "%.0f+-%.0f", mean(), sigma());
    return std::string(buf);
}

} // namespace pilon
