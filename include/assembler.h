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

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "pileup.h"

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
    std::string cigar;
};

class Assembler {
public:
    static int K;
    static constexpr int minDepth = 5;
    static constexpr int minExtend = 20;
    static constexpr int maxBranches = 5;
    static constexpr int minNovel = 200;
    static constexpr int minNovelPct = 50;

    using Kmer = std::string;
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

    std::string kmerPathString(const std::vector<std::string>& kmers, bool prependLength = false);
    void noteKmerLoop(int loopIndex, const std::vector<std::string>& kmers);

    std::vector<std::vector<std::string>>
    kmerPathsForward(std::vector<std::string> kmersIn, int branches = 0);
};

} // namespace pilon

#endif // PILON_ASSEMBLER_H
