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
#include <chrono>

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

bool BamFile::open(int ioThreads) {
    htsFile_ = hts_open(path_.c_str(), "rb");
    if (!htsFile_) {
        std::cerr << "Error: Cannot open BAM file: " << path_ << std::endl;
        return false;
    }

    if (ioThreads > 1) {
        hts_set_threads(htsFile_, ioThreads);
        hts_set_cache_size(htsFile_, Pilon::cacheSizeMb * 1024 * 1024);
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
    
    if (!index_) {
        std::cerr << "Error: BAM file " << path_ << " does not appear to have an index (.bai or .csi)" << std::endl;
        exit(1);
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

void BamFile::shareScanState(const BamFile& parent) {
    sharedStrayMap_ = &parent.strayMateMap_;      // pointer share, zero-copy
    insertStatsFR_ = parent.insertStatsFR_;
    insertStatsRF_ = parent.insertStatsRF_;
    insertStatsUnpaired_ = parent.insertStatsUnpaired_;
    bamType_ = parent.bamType_;
}

std::unordered_set<std::string> BamFile::getSeqNames() const {
    std::unordered_set<std::string> names;
    if (!header_) return names;

    int nseq = sam_hdr_nref(header_);
    for (int i = 0; i < nseq; i++) {
        const char* name = sam_hdr_tid2name(header_, i);
        if (name) {
            names.insert(name);
        }
    }
    return names;
}

std::vector<std::pair<std::string, int>> BamFile::getRefNamesAndLengths() const {
    std::vector<std::pair<std::string, int>> result;
    if (!header_) return result;
    int nseq = sam_hdr_nref(header_);
    for (int i = 0; i < nseq; i++) {
        const char* name = sam_hdr_tid2name(header_, i);
        uint32_t len = sam_hdr_tid2len(header_, i);
        if (name) result.emplace_back(name, static_cast<int>(len));
    }
    return result;
}

// Helper: check if mate is unmapped
static bool isMateUnmapped(const bam1_t* bam) {
    return bam->core.mtid < 0;
}

// Helper: convert a bam1_t to BamRead
static BamRead bamToRead(const bam1_t* bam, const bam_hdr_t* hdr) {
    BamRead read;
    read.readName = std::string(bam_get_qname(bam));
    
    // Convert encoded bases (4 bits per base) to ASCII
    std::string decoded;
    decoded.reserve(bam->core.l_qseq);
    for (int i = 0; i < bam->core.l_qseq; i++) {
        int base = bam_seqi(bam_get_seq(bam), i);
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

    // Store raw CIGAR ops (htslib uint32 encoding: len<<4 | op); avoid string
    // construction and later re-parsing (scaffold consumes these directly).
    uint32_t ncigar = bam->core.n_cigar;
    const uint32_t* cigarOps = bam_get_cigar(bam);
    read.cigar.assign(cigarOps, cigarOps + ncigar);

    return read;
}

// Helper: query BAM for reads in a region using htslib iterator
static std::vector<BamRead> queryRegion(const htsFile* fp, const bam_hdr_t* hdr,
                                         hts_idx_t* idx, const std::string& refName,
                                         int start, int stop) {
    std::vector<BamRead> reads;
    if (!fp || !hdr) return reads;

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

// Helper: add a base to pileup (matching Scala PileUpRegion.add)
// Scala: if (pair) { pileup.add + baseCount++ } else { badPair++ }
static void addPileup(GenomeRegion& region, int locus, char base, int qual, int mq, bool pair) {
    if (region.inRegion(locus)) {
        if (pair) {
            region.pileUpRegion(region.index(locus)).add(base, qual, mq);
            region.baseCount++;
        } else {
            region.pileUpRegion(region.index(locus)).badPair++;
        }
    }
}

// Helper: remove a base from pileup (matching Scala PileUpRegion.remove)
// Scala: if (pair) { no-op } else { badPair-- }
static void removePileup(GenomeRegion& region, int locus, char /*base*/, int /*qual*/, int /*mq*/, bool pair) {
    if (!region.inRegion(locus)) return;
    if (!pair) {
        region.pileUpRegion(region.index(locus)).badPair--;
    }
}

// Helper: physical coverage increment
static int physCovIncr(GenomeRegion& region, int aStart, int aEnd, int iSize, bool paired, bool valid) {
    if (!valid || (paired && iSize <= 0)) return 0;
    
    int start, end;
    if (!paired) {
        start = std::min(aStart, aEnd);
        end = std::max(aStart, aEnd);
    } else if (iSize > 0) {
        start = aStart;
        end = aStart + iSize;
    } else {
        // aEnd from bam_endpos is already 0-based exclusive, no +1 needed (matching Scala)
        end = aEnd;
        start = end + iSize; // iSize is negative
    }
    
    int insertSize = end - start;
    
    if (region.inRegion(start)) {
        region.pileUpRegion(region.index(start)).physCov++;
        region.pileUpRegion(region.index(start)).insertSize += insertSize;
    } else if (region.beforeRegion(start) && !region.beforeRegion(end)) {
        region.physCovStart++;
        region.insertSizeStart += insertSize;
    }
    
    if (region.inRegion(end)) {
        region.pileUpRegion(region.index(end)).physCov--;
        region.pileUpRegion(region.index(end)).insertSize -= insertSize;
    }
    
    return insertSize;
}

// Core: process a single read and add to pileups (matching Scala PileUpRegion.addRead)
static int addReadToPileup(GenomeRegion& region, const bam1_t* bam, int longReadType) {
    int length = bam->core.l_qseq;
    if (length == 0) return 0;

    const uint8_t* quals = bam_get_qual(bam);
    int mq = bam->core.qual;
    bool paired = (bam->core.flag & BAM_FPAIRED) != 0;
    
    // Valid read check
    bool valid = (mq >= Pilon::minMq) && 
                 ((!paired) || ((bam->core.flag & BAM_FPROPER_PAIR) && 
                   (bam->core.tid == bam->core.mtid)));
    
    int insert = bam->core.isize;
    int aStart = bam->core.pos;
    int aEnd = bam_endpos(bam);
    
    // Default qualities if none present (htslib stores absent quals as 0xFF or 0)
    bool noQuals = (quals[0] == 0 || quals[0] == 0xFF);
    if (noQuals && length > 0) {
        for (int i = 1; i < length; i++) {
            if (quals[i] != quals[0]) {
                noQuals = false;
                break;
            }
        }
    }
    // Reusable scratch buffers (thread-safe across worker threads via thread_local)
    static thread_local std::vector<uint8_t> qualsVec;
    qualsVec.resize(length);
    if (noQuals) {
        std::fill(qualsVec.begin(), qualsVec.end(), Pilon::defaultQual);
    } else {
        std::copy(quals, quals + length, qualsVec.begin());
    }
    
    // Trusted flank
    int trustedFlank = Pilon::flank;
    auto trusted = [trustedFlank, length](int offset) -> bool {
        return offset >= trustedFlank && (length - trustedFlank) > offset;
    };
    
    // Count clipped bases
    uint32_t ncigar = bam->core.n_cigar;
    const uint32_t* cigar = bam_get_cigar(bam);
    int clippedBases = 0;
    for (uint32_t i = 0; i < ncigar; i++) {
        if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
            clippedBases += bam_cigar_oplen(cigar[i]);
        }
    }
    
    // Adjust MQ proportionally to fraction of bases clipped
    int adjMq = Utils::roundDiv(mq * (length - clippedBases), length);
    int indelMq = longReadType > 0 ? std::min(adjMq, 8) : adjMq;
    
    // Decode bases to chars (htslib stores bases as 4-bit interleaved pairs)
    // Lookup table indexed by bam_seqi code: 1=A, 2=C, 4=G, 8=T, otherwise N
    static const char baseTable[16] = {'N','A','C','N','G','N','N','N','T','N','N','N','N','N','N','N'};
    static thread_local std::vector<char> baseChars;
    baseChars.resize(length);
    {
        const uint8_t* seq = bam_get_seq(bam);
        for (int i = 0; i < length; i++) {
            baseChars[i] = baseTable[bam_seqi(seq, i) & 15];
        }
    }
    
    // Parse CIGAR and add to pileups
    int readOffset = 0;
    int refOffset = 0;
    
    for (uint32_t i = 0; i < ncigar; i++) {
        int len = bam_cigar_oplen(cigar[i]);
        int op = bam_cigar_op(cigar[i]);
        int locus = aStart + refOffset;
        
        switch (op) {
            case BAM_CINS: { // Insertion
                if (valid && trusted(readOffset) && region.inRegion(locus)) {
                    // Get insertion sequence
                    std::vector<uint8_t> insertion(len);
                    for (int j = 0; j < len; j++) {
                        insertion[j] = static_cast<uint8_t>(baseChars[readOffset + j]);
                    }
                    
                    // Reposition insertion (matching Scala: while iloc > 1 && refBases(iloc - 2) == insertion(len - 1))
                    // Scala is 1-based, C++ is 0-based: iloc > 1 (1-based) = iloc > 0 (0-based)
                    // refBases(iloc - 2) in 1-based = refBase(iloc - 2) in 0-based
                    int iloc = locus;
                    while (iloc > 0 && region.inRegion(iloc - 1) &&
                           region.refBase(iloc - 1) == static_cast<char>(insertion.back())) {
                        iloc--;
                        // Rotate insertion
                        char last = insertion.back();
                        insertion.pop_back();
                        insertion.insert(insertion.begin(), last);
                    }
                    
                    // Skip if in homopolymer run (for long reads)
                    if (!(longReadType > 0 && region.homoRun(iloc - region.start) >= 4)) {
                        int idx = region.index(iloc);
                        if (idx >= 0 && idx < region.size()) {
                            uint8_t qual = qualsVec[readOffset];
                            region.pileUpRegion(idx).addInsertion(insertion, qual, indelMq);
                        }
                    }
                }
                break;
            }
            case BAM_CDEL: { // Deletion
                if (valid && trusted(readOffset) && region.inRegion(locus) && 
                    region.inRegion(locus + len - 1)) {
                    int dloc = locus;
                    int rloc = readOffset;
                    
                    // Slide deletion (matching Scala: while dloc > 1 && rloc > 0 && refBases(dloc - 2) == refBases(dloc + len - 2))
                    // Scala: dloc is 1-based, refBases is 0-based
                    // C++: dloc is 0-based, refBase() is 0-based
                    // dloc > 1 (1-based) = dloc > 0 (0-based)
                    // refBases(dloc - 2) (1-based dloc) = refBase(dloc - 1) (0-based dloc)
                    // refBases(dloc + len - 2) (1-based dloc) = refBase(dloc + len - 1) (0-based dloc)
                    while (dloc > 0 && rloc > 0 && 
                           region.inRegion(dloc - 1) && region.inRegion(dloc + len - 1) &&
                           region.refBase(dloc - 1) == region.refBase(dloc + len - 1)) {
                        dloc--;
                        rloc--;
                        char base = baseChars[rloc];
                        uint8_t qual = qualsVec[rloc];
                        
                        if (trusted(rloc) && region.inRegion(dloc)) {
                            removePileup(region, dloc, base, qual, adjMq, valid);
                            if (region.inRegion(dloc + len)) {
                                addPileup(region, dloc + len, base, qual, adjMq, valid);
                            }
                        }
                    }
                    
                    // Skip if in homopolymer run (for long reads)
                    // Scala: also check nanoporeExclude for nanopore reads
                    bool skipDeletion = (longReadType > 0 && region.homoRun(dloc - region.start) >= 4);
                    if (!skipDeletion && longReadType == BamFile::nanoporeLongRead) {
                        skipDeletion = region.nanoporeExclude(dloc - region.start);
                    }
                    if (!skipDeletion) {
                        // Get deletion sequence from reference
                        // Scala: refBases.slice(dloc - 1, dloc + len - 1) where dloc is 1-based
                        // C++: dloc is 0-based, so range is [dloc, dloc + len) = len elements
                        std::vector<uint8_t> deletion;
                        for (int j = dloc; j < dloc + len; j++) {
                            if (region.inRegion(j)) {
                                deletion.push_back(static_cast<uint8_t>(region.refBase(j)));
                            }
                        }
                        
                        if (!deletion.empty()) {
                            int idx = region.index(dloc);
                            if (idx >= 0 && idx < region.size()) {
                                uint8_t qual = qualsVec[readOffset];
                                region.pileUpRegion(idx).addDeletion(deletion, qual, indelMq);
                            }
                        }
                    }
                }
                break;
            }
            case BAM_CMATCH: // M, EQ, X
            case BAM_CEQUAL:
            case BAM_CDIFF: {
                for (int j = 0; j < len; j++) {
                    int rOff = readOffset + j;
                    if (trusted(rOff)) {
                        int locusPlus = locus + j;
                        char base = baseChars[rOff];
                        uint8_t qual = qualsVec[rOff];
                        // Scala: Nanopore CCGG motif exclude (set qual=0)
                        if (longReadType == BamFile::nanoporeLongRead &&
                            region.inRegion(locusPlus) &&
                            region.nanoporeExclude(region.index(locusPlus))) {
                            qual = 0;
                        }
                        addPileup(region, locusPlus, base, qual, adjMq, valid);
                    }
                }
                break;
            }
            case BAM_CSOFT_CLIP: { // Soft clip
                int clipStart = (readOffset == 0) ? locus - len : locus;
                int clipEnd = clipStart + len - 1;
                
                if (region.inRegion(clipStart)) {
                    region.pileUpRegion(region.index(clipStart)).clips++;
                }
                if (region.inRegion(clipEnd)) {
                    region.pileUpRegion(region.index(clipEnd)).clips++;
                }
                
                for (int j = 0; j < len; j++) {
                    int rOff = readOffset + j;
                    int locusPlus = clipStart + j;
                    char base = baseChars[rOff];
                    uint8_t qual = qualsVec[rOff];
                    if (region.inRegion(locusPlus)) {
                        addPileup(region, locusPlus, base, qual, adjMq, false);
                    }
                }
                break;
            }
            case BAM_CHARD_CLIP: // Hard clip - nothing to do
                break;
            case BAM_CREF_SKIP: // N - skip reference bases
                break;
            default:
                if (Pilon::verbose) {
                    std::cerr << "Unknown CIGAR op=" << bam_cigar_opchr(op) 
                              << " in read " << bam_get_qname(bam) << std::endl;
                }
                break;
        }
        
        // Update offsets - manually check CIGAR op types for older htslib
        bool consumesRead = (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CSOFT_CLIP || 
                             op == BAM_CEQUAL || op == BAM_CDIFF);
        bool consumesRef = (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || 
                            op == BAM_CEQUAL || op == BAM_CDIFF);
        if (consumesRead) readOffset += len;
        if (consumesRef) refOffset += len;
    }
    
    region.readCount++;
    int insertSz = physCovIncr(region, aStart, aEnd, insert, paired, valid);
    return insertSz;
}

double BamFile::process(GenomeRegion& region, int printInterval) {
    if (!htsFile_) {
        if (!open()) return 0;
    }

    int tid = sam_hdr_name2tid(header_, region.name.c_str());
    if (tid < 0) return 0;

    hts_itr_t* iter = nullptr;
    if (index_) {
        int contigLen = sam_hdr_tid2len(header_, tid);
        int queryStop = std::min(region.stop + 10000, contigLen);
        iter = sam_itr_queryi(index_, tid, std::max(0, region.start - 10000), queryStop);
    } else {
        int contigLen = sam_hdr_tid2len(header_, tid);
        int queryStop = std::min(region.stop + 10000, contigLen);
        iter = sam_itr_queryi(nullptr, tid, std::max(0, region.start - 10000), queryStop);
    }
    if (!iter) return 0;

    bam1_t* bam = bam_init1();
    int lastLoc = 0;
    long long processedReads = 0;
    auto pileStart = std::chrono::steady_clock::now();
    auto pileLast = pileStart;

    int ret;
    while ((ret = sam_itr_next(htsFile_, iter, bam)) > 0) {
        int loc = bam->core.pos;
        if (Pilon::verbose && printInterval > 0 && loc > lastLoc + printInterval) {
            lastLoc = printInterval * (loc / printInterval);
            std::cout << "..." << lastLoc << std::flush;
        }

        if (validateRead(bam)) {
            int computedInsert = addReadToPileup(region, bam, longReadType_);
            bool rc = (bam->core.flag & BAM_FREVERSE) != 0;
            addInsert(computedInsert, rc);
        }

        processedReads++;
        if (processedReads % 100000 == 0) {
            auto now = std::chrono::steady_clock::now();
            double recent = std::chrono::duration<double>(now - pileLast).count();
            if (recent >= 2.0 && !Pilon::verbose) {
                double elapsed = std::chrono::duration<double>(now - pileStart).count();
                double readsPerSec = processedReads / elapsed;
                printf("\r  pileup %s:%d-%d... %lld reads | %.1fM/s | elapsed %.1fs  ",
                       region.name.c_str(), region.start, region.stop,
                       (long long)processedReads, readsPerSec / 1e6, elapsed);
                fflush(stdout);
                pileLast = now;
            }
        }
    }

    // Clear progress line for non-verbose mode
    if (!Pilon::verbose) {
        auto pileEnd = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(pileEnd - pileStart).count();
        printf("\r  pileup %s:%d-%d... %lld reads in %.1fs\n",
               region.name.c_str(), region.start, region.stop,
               (long long)processedReads, elapsed);
    }

    if (ret < -1) {
        std::cerr << "Error reading BAM file" << std::endl;
    }

    hts_itr_destroy(iter);
    bam_destroy1(bam);

    // Compute mean coverage for this region (matching Scala return value)
    double meanCoverage = 0.0;
    if (region.size() > 0) {
        long long totalCov = 0;
        for (int i = 0; i < region.size(); i++)
            totalCov += region.pileUpRegion(i).depth();
        meanCoverage = static_cast<double>(totalCov) / region.size();
    }
    return meanCoverage;
}

long long BamFile::estimateTotalReads() const {
    if (!index_ || !header_) return 0;
    long long total = 0;
    int nref = sam_hdr_nref(header_);
    for (int tid = 0; tid < nref; tid++) {
        uint64_t mapped = 0, unmapped = 0;
        if (hts_idx_get_stat(index_, tid, &mapped, &unmapped) >= 0) {
            total += static_cast<long long>(mapped);
        }
    }
    if (total > 0) return total;
    uint64_t noco = hts_idx_get_n_no_coor(index_);
    return static_cast<long long>(noco);
}

void BamFile::scan(const std::unordered_set<std::string>& seqsOfInterest) {
    if (!htsFile_) {
        if (!open()) return;
    }

    long long estTotal = estimateTotalReads();

    auto startTime = std::chrono::steady_clock::now();
    long long lastReport = 0;
    auto lastTime = startTime;

    bam1_t* bam = bam_init1();
    int ret;
    long long totalReads = 0;

    while ((ret = sam_read1(htsFile_, header_, bam)) > 0) {
        totalReads++;
        if (!validateRead(bam)) {
            filtered_++;
        } else if (bam->core.flag & BAM_FUNMAP) {
            unmapped_++;
        } else {
            mapped_++;
            if (bam->core.flag & BAM_FPAIRED) {
                if (bam->core.flag & BAM_FPROPER_PAIR) {
                    proper_++;
                    addInsert(bam->core.isize,
                             (bam->core.flag & BAM_FREVERSE) != 0);
                } else if (Pilon::strays && !isMateUnmapped(bam)) {
                    const char* refName = sam_hdr_tid2name(header_, bam->core.tid);
                    const char* mateRef = sam_hdr_tid2name(header_, bam->core.mtid);
                    bool refInInterest = refName && seqsOfInterest.count(refName);
                    bool mateInInterest = mateRef && seqsOfInterest.count(mateRef);
                    if (refInInterest || mateInInterest) {
                        BamRead read = bamToRead(bam, header_);
                        strayMateMap_.addRead(read);
                    }
                }
            } else {
                addInsert(bam->core.l_qseq,
                         (bam->core.flag & BAM_FREVERSE) != 0, true);
            }
        }

        if (totalReads - lastReport >= 500000) {
            auto now = std::chrono::steady_clock::now();
            double recent = std::chrono::duration<double>(now - lastTime).count();
            if (recent >= 2.0) {
                double elapsed = std::chrono::duration<double>(now - startTime).count();
                double readsPerSec = totalReads / elapsed;
                if (estTotal > 0) {
                    int pct = static_cast<int>(totalReads * 100 / estTotal);
                    double eta = (estTotal - totalReads) / readsPerSec;
                    printf("\r  Scanning %s... %d%% %lld/%lld reads | %.1fM/s | ETA %.1fs  ",
                           path_.c_str(), pct, (long long)totalReads, (long long)estTotal,
                           readsPerSec / 1e6, eta);
                } else {
                    printf("\r  Scanning %s... %lld reads | %.1fM/s | elapsed %.1fs  ",
                           path_.c_str(), (long long)totalReads, readsPerSec / 1e6, elapsed);
                }
                fflush(stdout);
                lastReport = totalReads;
                lastTime = now;
            }
        }
    }

    bam_destroy1(bam);

    auto endTime = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(endTime - startTime).count();
    printf("\r  Scanning %s... done (%lld reads in %.1fs)\n",
           path_.c_str(), (long long)totalReads, elapsed);

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
    auto reads = queryRegion(htsFile_, header_, index_, region.name, region.start, region.stop);
    // Filter reads matching Scala validateRead (QC-fail, duplicate, secondary alignment)
    reads.erase(std::remove_if(reads.begin(), reads.end(),
        [](const BamRead& r) {
            return (!Pilon::nonPf && r.failsVendorQC) ||
                   (!Pilon::duplicates && r.duplicate) ||
                   r.secondary;
        }), reads.end());
    return reads;
}

std::vector<BamRead> BamFile::recruitFlankReads(const Region& region) const {
    Region flanks = flankRegion(region);
    std::vector<BamRead> reads = readsInRegion(flanks);
    
    // Handle stray mates if Pilon::strays is set (matching Scala: if (bamType != "unpaired"))
    if (Pilon::strays && bamType_ != "unpaired") {
        std::unordered_map<std::string, BamRead> readMap1;
        std::unordered_map<std::string, BamRead> readMap2;
        for (const auto& r : reads) {
            if (r.firstOfPair) {
                readMap1[r.readName] = r;
            } else {
                readMap2[r.readName] = r;
            }
        }
        
        for (const auto& [name, read] : readMap1) {
            if (readMap2.find(name) == readMap2.end()) {
                BamRead* mate = getStrayMap().lookupMate(read);
                if (mate) reads.push_back(*mate);
            }
        }
        for (const auto& [name, read] : readMap2) {
            if (readMap1.find(name) == readMap1.end()) {
                BamRead* mate = getStrayMap().lookupMate(read);
                if (mate) reads.push_back(*mate);
            }
        }
    }
    
    return reads;
}

std::vector<BamRead> BamFile::recruitBadMates(const Region& region) const {
    int midpoint = (region.start + region.stop) / 2;
    Region flanks = flankRegion(region);
    std::vector<BamRead> reads = readsInRegion(flanks);

    // Build local mate maps (matching Scala MateMap)
    std::unordered_map<std::string, const BamRead*> readMap1;
    std::unordered_map<std::string, const BamRead*> readMap2;
    for (const auto& r : reads) {
        if (r.firstOfPair) readMap1[r.readName] = &r;
        else readMap2[r.readName] = &r;
    }

    // Add stray mates (matching Scala mateMap.addStrays)
    if (Pilon::strays) {
        for (const auto& [name, r] : readMap1) {
            if (readMap2.find(name) == readMap2.end()) {
                BamRead* mate = getStrayMap().lookupMate(*r);
                if (mate) readMap2[name] = mate;
            }
        }
        for (const auto& [name, r] : readMap2) {
            if (readMap1.find(name) == readMap1.end()) {
                BamRead* mate = getStrayMap().lookupMate(*r);
                if (mate) readMap1[name] = mate;
            }
        }
    }

    int frPctVal = pctFR();

    // Filter for anchored reads with good orientation (matching Scala)
    std::vector<BamRead> mates;
    for (const auto& [name, r1] : readMap1) {
        auto it = readMap2.find(name);
        if (it == readMap2.end()) continue;
        const auto& r2 = *it->second;

        // Check r1 is not unmapped and not a proper pair
        if (r1->unmapped || r1->properPair) continue;
        
        const auto& r = *r1;  // the anchored read
        bool rc = r.negativeStrand;
        int start = r.alignmentStart;
        int end = r.alignmentEnd;
        bool before = start < midpoint;
        bool after = end > midpoint;

        // FR orientation: anchored before midpoint on + strand, or after on - strand
        bool frOrientation = (before && !rc) || (after && rc);

        bool goodOrientation;
        if (frPctVal > 100 - minOrientationPct) goodOrientation = frOrientation;
        else if (frPctVal < minOrientationPct) goodOrientation = !frOrientation;
        else goodOrientation = true;

        if (goodOrientation && ((flanks.start <= start && start < flanks.stop) ||
                                 (flanks.start <= end && end < flanks.stop)))
            mates.push_back(r2);
    }

    // Reverse pass: second-of-pair as anchor (matching Scala pairs() which emits both orders)
    for (const auto& [name, r2] : readMap2) {
        auto it = readMap1.find(name);
        if (it == readMap1.end()) continue;
        const auto& r1 = *it->second;

        if (r2->unmapped || r2->properPair) continue;
        
        const auto& r = *r2;  // the anchored read (second of pair)
        bool rc = r.negativeStrand;
        int start = r.alignmentStart;
        int end = r.alignmentEnd;
        bool before = start < midpoint;
        bool after = end > midpoint;

        bool frOrientation = (before && !rc) || (after && rc);

        bool goodOrientation;
        if (frPctVal > 100 - minOrientationPct) goodOrientation = frOrientation;
        else if (frPctVal < minOrientationPct) goodOrientation = !frOrientation;
        else goodOrientation = true;

        if (goodOrientation && ((flanks.start <= start && start < flanks.stop) ||
                                 (flanks.start <= end && end < flanks.stop)))
            mates.push_back(r1);  // the other read (first of pair)
    }

    return mates;
}

// Get unmapped reads for novel contig assembly
std::vector<BamRead> BamFile::getUnalignedReads() const {
    std::vector<BamRead> results;
    htsFile* fp = hts_open(path_.c_str(), "r");
    if (!fp) return results;
    bam_hdr_t* hdr = sam_hdr_read(fp);
    if (!hdr) { hts_close(fp); return results; }

    bam1_t* bam = bam_init1();
    int ret;
    while ((ret = sam_read1(fp, hdr, bam)) > 0) {
        if (!Pilon::nonPf && (bam->core.flag & BAM_FQCFAIL)) continue;
        if (!Pilon::duplicates && (bam->core.flag & BAM_FDUP)) continue;
        if (bam->core.flag & BAM_FSECONDARY) continue;
        if (bam->core.flag & BAM_FUNMAP) {
            results.push_back(bamToRead(bam, hdr));
        }
    }
    bam_destroy1(bam);
    bam_hdr_destroy(hdr);
    hts_close(fp);
    return results;
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
    const int huge = 5 * getMaxInsertSizes()[bamType_];
    
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
    return Utils::pct(insertStatsFR_.count, total);
}

int BamFile::pctRF() const {
    long long total = insertStatsFR_.count + insertStatsRF_.count + insertStatsUnpaired_.count;
    if (total == 0) return 0;
    return Utils::pct(insertStatsRF_.count, total);
}

int BamFile::pctUnpaired() const {
    long long total = insertStatsFR_.count + insertStatsRF_.count + insertStatsUnpaired_.count;
    if (total == 0) return 0;
    return Utils::pct(insertStatsUnpaired_.count, total);
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
    if (count <= 1000) return 0;
    return static_cast<int>(std::round(mean() + 3.0 * sigma()));
}

std::string BamFile::InsertSizeStats::toString() const {
    char buf[64];
    snprintf(buf, sizeof(buf), "%.0f+-%.0f", mean(), sigma());
    return std::string(buf);
}

} // namespace pilon
