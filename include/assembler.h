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

#ifndef PILON_ASSEMBLER_H
#define PILON_ASSEMBLER_H

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "pileup.h"

namespace pilon {

// Packed k-mer representation for the assembler.
//
// Most reads/sequences are pure ACGT, encoded as 2 bits per base in a 94-bit
// integer (K=47). The rare k-mers containing N (or any non-ACGT base) fall back
// to their literal string and are tagged with hasN so they remain distinct keys,
// exactly matching the Scala string-keyed behavior. Non-ACGT kmers never reach
// minDepth and never enter the graph, so the fallback is essentially inert.
struct Kmer {
    uint64_t hi = 0;   // bits 64..(2*K-1)
    uint64_t lo = 0;   // bits 0..63
    bool hasN = false; // true => this k-mer is represented by `s`
    std::string s;     // used only when hasN is set

    bool operator==(const Kmer& o) const {
        if (hasN != o.hasN) return false;
        if (hasN) return s == o.s;
        return hi == o.hi && lo == o.lo;
    }
    bool operator!=(const Kmer& o) const { return !(*this == o); }
};

} // namespace pilon

namespace std {
template<> struct hash<pilon::Kmer> {
    size_t operator()(const pilon::Kmer& k) const noexcept {
        if (k.hasN) return std::hash<std::string>()(k.s);
        uint64_t x = k.hi ^ (k.lo + 0x9e3779b97f4a7c15ULL + (k.hi << 6) + (k.hi >> 2));
        x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL; x ^= x >> 27; x *= 0x94d049bb133111ebULL; x ^= x >> 31;
        return static_cast<size_t>(x);
    }
};
}

namespace pilon {

struct BamRead {
    std::string readName;
    std::string bases;
    std::vector<uint8_t> quals;
    int mappingQuality;
    int alignmentStart;
    int alignmentEnd;
    bool negativeStrand;
    bool paired;
    bool properPair;
    bool firstOfPair;
    bool unmapped;
    bool mateUnmapped;
    bool duplicate;
    bool secondary;
    bool failsVendorQC;
    std::string referenceName;
    std::string mateReferenceName;
    int mateAlignmentStart;
    int inferredInsertSize;
    std::vector<uint32_t> cigar;  // raw htslib CIGAR ops (len<<4 | op), avoids re-parsing
};

class Assembler {
public:
    static int K;
    static constexpr int minDepth = 5;
    static constexpr int minExtend = 20;
    static constexpr int maxBranches = 5;
    static constexpr int minNovel = 200;
    static constexpr int minNovelPct = 50;

    using KmerPileup = std::unordered_map<Kmer, PileUp>;
    using KmerGraph = std::unordered_map<Kmer, Kmer>;

    Assembler(int minDepth_ = minDepth);

    void addReads(const std::vector<BamRead>& reads);
    void addPair(const BamRead& r1, const BamRead& r2);
    void addRead(const BamRead& r);
    void addSeq(const std::string& bases);
    void addGraphSeq(const std::string& bases);
    void addGraphSeqs(const std::vector<std::string>& seqs);

    void buildGraph();
    void prunePileups(int minCount = minDepth);

    std::vector<std::string> pathsForward(const std::string& startingKmer);
    std::vector<std::string> pathsReverse(const std::string& startingKmer);
    std::vector<std::string> tryForward(const std::string& anchor);
    std::vector<std::string> tryReverse(const std::string& anchor);

    // Returns {forwardPaths, reversePaths, loopSequence}
    std::tuple<std::vector<std::string>, std::vector<std::string>, std::string>
    multiBridge(const std::string& left, const std::string& right);

    std::vector<std::string> novel(Assembler& ref);

    // Accessors
    const KmerGraph& getKGraph() const { return kGraph; }
    const KmerPileup& getPileups() const { return pileups; }

    long long nReads() const { return nReads_; }
    long long nBases() const { return nBases_; }

    std::string toString() const;

private:
    int minDepth_;
    KmerPileup pileups;
    KmerGraph kGraph;
    KmerGraph altGraph;
    long long nReads_;
    long long nBases_;
    int loopLength_;
    std::string loopSequence_;

    void addToPileups(const std::string& bases, const std::vector<uint8_t>& quals, int mq);
    void graphSeq(const std::string& bases);
    static void addLink(KmerGraph& g, const Kmer& k1, const Kmer& k2, int weight);

    std::string kmerPathString(const std::vector<Kmer>& kmers, bool prependLength = false);
    void noteKmerLoop(int loopIndex, const std::vector<Kmer>& kmers);

    std::vector<std::vector<Kmer>>
    kmerPathsForward(std::vector<Kmer> kmersIn, int branches = 0);
};

} // namespace pilon

#endif // PILON_ASSEMBLER_H
