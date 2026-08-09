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

#include "scaffold.h"
#include "pilon.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <set>
#include <map>

namespace pilon {

// =============================================================================
// MatePair (matching Scala class MatePair)
// =============================================================================

MatePair::MatePair(const BamRead& r1, const BamRead& r2,
                   const std::unordered_map<std::string, int>& refNameToIdx) {
    // Scala: scaffold(read) = read.getReferenceIndex
    auto scaffold = [&](const BamRead& r) -> int {
        auto it = refNameToIdx.find(r.referenceName);
        return it != refNameToIdx.end() ? it->second : -1;
    };

    // Scala: coord(read)
    //   fw coords are negative, rc are positive
    //   => if negativeStrand: alignmentEnd, else: -alignmentStart
    auto coord = [](const BamRead& r) -> int {
        if (r.negativeStrand) return r.alignmentEnd;
        else return -r.alignmentStart;
    };

    int s1 = scaffold(r1);
    int s2 = scaffold(r2);
    int c1 = coord(r1);
    int c2 = coord(r2);

    // Order lowest first
    if ((s1 < s2) || (s1 == s2 && c1 < c2)) {
        scaffold1 = s1; coord1 = c1;
        scaffold2 = s2; coord2 = c2;
    } else {
        scaffold1 = s2; coord1 = c2;
        scaffold2 = s1; coord2 = c1;
    }

    mq = std::min(r1.mappingQuality, r2.mappingQuality);
}

int64_t MatePair::longCoord1() const {
    // Scala: longCoord(s: Int, c: Int): Long = (s.toLong << 32) | (c.toLong & mask32bit)
    return (static_cast<int64_t>(scaffold1) << 32) | (static_cast<int64_t>(coord1) & 0xFFFFFFFFL);
}

int64_t MatePair::longCoord2() const {
    return (static_cast<int64_t>(scaffold2) << 32) | (static_cast<int64_t>(coord2) & 0xFFFFFFFFL);
}

int64_t MatePair::distance() const {
    // Scala: if (sameScaffold) (coord1 + coord2)
    //        else ((scaffold1 & mask16bit).toLong << 48) | ((scaffold2 & mask16bit) << 32) | ((coord1 + coord2) & mask32bit)
    int64_t sum = static_cast<int64_t>(coord1) + static_cast<int64_t>(coord2);
    if (sameScaffold()) {
        return sum;
    } else {
        int64_t s1 = static_cast<int64_t>(scaffold1 & 0xFFFF);
        int64_t s2 = static_cast<int64_t>(scaffold2 & 0xFFFF);
        return (s1 << 48) | (s2 << 32) | (sum & 0xFFFFFFFFL);
    }
}

bool MatePair::isRunt(int maxInsert) const {
    // Scala: if (sameScaffold) { val d = distance; (d < maxInsert) && (d > -maxInsert) } else false
    if (!sameScaffold()) return false;
    int64_t d = distance();
    return (d < maxInsert) && (d > -maxInsert);
}

// =============================================================================
// LinkCluster (matching Scala class LinkCluster)
// =============================================================================

LinkCluster::LinkCluster(const std::vector<MatePair>& mates,
                         const std::vector<std::string>& refNames,
                         const std::vector<int>& refLengths,
                         int sigma)
    : nLinks(static_cast<int>(mates.size())), sigma_(sigma) {
    scaffold1 = mates[0].scaffold1;
    minCoord1 = mates[0].coord1;
    maxCoord1 = mates[0].coord1;
    scaffold2 = mates[0].scaffold2;
    minCoord2 = mates[0].coord2;
    maxCoord2 = mates[0].coord2;
    int mqSum = mates[0].mq;
    for (size_t i = 1; i < mates.size(); i++) {
        minCoord1 = std::min(minCoord1, mates[i].coord1);
        maxCoord1 = std::max(maxCoord1, mates[i].coord1);
        minCoord2 = std::min(minCoord2, mates[i].coord2);
        maxCoord2 = std::max(maxCoord2, mates[i].coord2);
        mqSum += mates[i].mq;
    }
    mq = mqSum / static_cast<int>(mates.size());
    seqName1 = (scaffold1 >= 0 && scaffold1 < static_cast<int>(refNames.size())) ? refNames[scaffold1] : "?";
    seqName2 = (scaffold2 >= 0 && scaffold2 < static_cast<int>(refNames.size())) ? refNames[scaffold2] : "?";
    size1 = (scaffold1 >= 0 && scaffold1 < static_cast<int>(refLengths.size())) ? refLengths[scaffold1] : 0;
    size2 = (scaffold2 >= 0 && scaffold2 < static_cast<int>(refLengths.size())) ? refLengths[scaffold2] : 0;
}

std::string LinkCluster::reportCoord1() const {
    // Scala: val dir = if (coord < 0) "+" else "-"
    //        "%s:%d%s".format(name, coord.abs, dir)
    std::string dir = (minCoord1 < 0) ? "+" : "-";
    return seqName1 + ":" + std::to_string(std::abs(minCoord1)) + dir;
}

std::string LinkCluster::reportCoord2() const {
    std::string dir = (minCoord2 < 0) ? "+" : "-";
    return seqName2 + ":" + std::to_string(std::abs(minCoord2)) + dir;
}

// =============================================================================
// EndAlignment (matching Scala class EndAlignment)
// =============================================================================

EndAlignment::EndAlignment(const BamRead& read, const std::string& cName, int cLen)
    : readName(read.readName), contigName(cName), contigLength(cLen),
      alignStart(read.alignmentStart), alignEnd(read.alignmentEnd),
      rc(read.negativeStrand) {
    // Compute unclipped coordinates from CIGAR soft clips
    int us = read.alignmentStart;
    int ue = read.alignmentEnd;
    const std::vector<uint32_t>& cig = read.cigar;
    if (!cig.empty()) {
        // First CIGAR op: left soft clip
        if (bam_cigar_op(cig.front()) == BAM_CSOFT_CLIP) {
            us -= static_cast<int>(bam_cigar_oplen(cig.front()));
        }
        // Last CIGAR op: right soft clip
        if (bam_cigar_op(cig.back()) == BAM_CSOFT_CLIP) {
            ue += static_cast<int>(bam_cigar_oplen(cig.back()));
        }
    }
    unclippedStart = us;
    unclippedEnd = ue;
    // Scala uses: align.getUnclippedStart, align.getUnclippedEnd
    // These come from SAMRecord which parses CIGAR soft clips
    // In our BamRead, we'd need to parse the cigar string
    // For now, approximate with alignmentStart/End
    int endLength = std::min(EndMaxLength, contigLength / 2);

    cantileveredLeft = unclippedStart < endLength &&
                       (unclippedEnd == alignEnd || true) &&
                       alignEnd > EndMinOverlap &&
                       (unclippedStart < 1 || unclippedStart == alignStart);
    cantileveredRight = unclippedEnd > contigLength - endLength &&
                        (unclippedStart == alignStart || true) &&
                        alignStart < contigLength - EndMinOverlap &&
                        (unclippedEnd > contigLength || unclippedEnd == alignEnd);
}

std::string EndAlignment::contigEndName() const {
    std::string str = contigName;
    if (cantileveredLeft) str += "L";
    if (cantileveredRight) str += "R";
    return str;
}

// =============================================================================
// EndAlignmentPair (matching Scala class EndAlignmentPair)
// =============================================================================

EndAlignmentPair::EndAlignmentPair(const EndAlignment& a1, const EndAlignment& a2)
    : ca1(a1), ca2(a2) {
    std::string e1 = ca1.contigEndName();
    std::string e2 = ca2.contigEndName();
    if (e1 >= e2) {
        ca1 = a1;
        ca2 = a2;
    } else {
        ca1 = a2;
        ca2 = a1;
    }
}

int EndAlignmentPair::impliedLength() const {
    // Scala: (Math.abs(ca2.unclippedStart - ca1.unclippedStart) +
    //          Math.abs(ca2.unclippedEnd - ca1.unclippedEnd) + 1) / 2
    return (std::abs(ca2.unclippedStart - ca1.unclippedStart) +
            std::abs(ca2.unclippedEnd - ca1.unclippedEnd) + 1) / 2;
}

// =============================================================================
// Scaffold namespace functions
// =============================================================================

void Scaffold::dumpCoords(const std::vector<MatePair>& coords) {
    if (!Pilon::debug) return;
    for (const auto& mp : coords) {
        // Print scaffold1:coord1 scaffold2:coord2 mq distance
        std::cout << mp.scaffold1 << ":" << mp.coord1 << " "
                  << mp.scaffold2 << ":" << mp.coord2 << " "
                  << mp.mq << " " << mp.distance() << std::endl;
    }
}

std::vector<std::vector<MatePair>> Scaffold::findClustersInternal(
    const std::vector<MatePair>& unsortedCoords, int window, int minCluster,
    std::function<int64_t(const MatePair&)> mpFunc) {

    auto coords = unsortedCoords; // copy
    std::sort(coords.begin(), coords.end(),
              [&mpFunc](const MatePair& a, const MatePair& b) {
                  return mpFunc(a) < mpFunc(b);
              });

    if (Pilon::debug) dumpCoords(coords);

    std::pair<int, int> best = {0, 0};
    std::vector<std::pair<int, int>> clusterRanges;

    for (int tail = 0; tail < static_cast<int>(coords.size()); tail++) {
        int64_t windowLimit = mpFunc(coords[tail]) + window;
        int head = tail + 1;
        while (head < static_cast<int>(coords.size()) && mpFunc(coords[head]) < windowLimit)
            head++;

        if (head - tail > minCluster) {
            if (head - tail > best.second - best.first)
                best = {tail, head};
        } else {
            if (best.second > 0) clusterRanges.push_back(best);
            best = {0, 0};
        }
    }
    if (best.second > 0) clusterRanges.push_back(best);

    std::vector<std::vector<MatePair>> result;
    for (const auto& range : clusterRanges) {
        result.emplace_back(coords.begin() + range.first, coords.begin() + range.second);
    }
    return result;
}

std::vector<LinkCluster> Scaffold::findClusters(
    std::vector<MatePair>& coords,
    const std::vector<std::string>& refNames,
    const std::vector<int>& refLengths,
    int sigma, int minCluster) {

    int windowSize = 4 * sigma;
    auto clustersByDistance = findClustersInternal(coords, windowSize, minCluster,
        [](const MatePair& mp) { return mp.distance(); });

    // Second pass: sub-cluster each distance cluster by longCoord1
    std::vector<LinkCluster> linkClusters;
    for (auto& cluster : clustersByDistance) {
        auto subClusters = findClustersInternal(cluster, windowSize, minCluster,
            [](const MatePair& mp) { return mp.longCoord1(); });
        for (auto& sc : subClusters) {
            linkClusters.emplace_back(sc, refNames, refLengths, sigma);
        }
    }

    std::sort(linkClusters.begin(), linkClusters.end(),
              [](const LinkCluster& a, const LinkCluster& b) {
                  return a.nLinks > b.nLinks;
              });

    if (Pilon::debug) {
        for (const auto& lc : linkClusters) {
            std::cout << "<LinkCluster " << lc.nLinks << " "
                      << lc.scaffold1 << ":" << lc.minCoord1 << "+" << lc.maxCoord1 << " "
                      << lc.scaffold2 << ":" << lc.minCoord2 << "+" << lc.maxCoord2 << " "
                      << lc.mq;
            if (lc.valid()) std::cout << " valid";
            if (lc.circular()) std::cout << " circular";
            std::cout << ">" << std::endl;
        }
    }
    return linkClusters;
}

void Scaffold::analyzeStrays(BamFile& bam) {
    auto pairs = bam.getStrayMap().pairs();
    int sigma = static_cast<int>(bam.insertSizeSigma());
    auto refInfo = bam.getRefNamesAndLengths();

    std::vector<std::string> refNames;
    std::vector<int> refLengths;
    int genomeSize = 0;
    for (const auto& ri : refInfo) {
        refNames.push_back(ri.first);
        refLengths.push_back(ri.second);
        genomeSize += ri.second;
    }

    // Build name-to-index map
    std::unordered_map<std::string, int> refNameToIdx;
    for (size_t i = 0; i < refNames.size(); i++)
        refNameToIdx[refNames[i]] = static_cast<int>(i);

    std::cout << "Analyzing large-scale structure using " << bam.path() << std::endl;
    if (Pilon::verbose) {
        std::cout << "analyzing strays in " << bam.path() << std::endl;
        std::cout << "genome size " << genomeSize << std::endl;
        std::cout << "max insert " << bam.maxInsertSize() << std::endl;
    }

    int innies = 0;
    int intra = 0;
    int ambig = 0;
    std::vector<MatePair> links;

    for (const auto& [r1, r2] : pairs) {
        if (!r2.firstOfPair) {  // r2 is second of pair
            MatePair mp(r1, r2, refNameToIdx);
            if (mp.isRunt(10000)) {
                innies++;
            } else if (mp.ambiguousPlacement()) {
                ambig++;
            } else {
                if (mp.sameScaffold()) intra++;
                links.push_back(mp);
            }
        }
    }

    int nLinks = static_cast<int>(links.size());
    int backgroundRate = static_cast<int>(nLinks * 4.0f * sigma / genomeSize);
    if (Pilon::debug) {
        std::cout << "ambig=" << ambig << " same=" << intra
                  << " innies=" << innies << " links=" << nLinks
                  << " background=" << backgroundRate << std::endl;
    }

    // minCluster: use 25 like Scala (TODO: backgroundRate)
    int minClusterSize = std::max(25, backgroundRate);
    auto clusters = findClusters(links, refNames, refLengths, sigma, minClusterSize);

    for (const auto& c : clusters) {
        if (!c.valid()) continue;
        auto c1 = c.reportCoord1();
        auto c2 = c.reportCoord2();
        if (c.circular()) {
            std::cout << "Circular element " << c.seqName1;
        } else if (c.scaffoldLink()) {
            std::cout << "Candidate scaffold link " << c1 << " to " << c2;
        } else if (c.rearrangement()) {
            std::cout << "Candidate rearrangement " << c1 << " connects to " << c2;
        } else {
            std::cout << "Other cluster " << c1 << " - " << c2;
        }
        if (c.sameOrientation()) std::cout << " reversed";
        std::cout << " (" << c.nLinks << " supporting links)" << std::endl;
    }
}

void Scaffold::analyze(std::vector<BamFile*>& bamFiles) {
    std::cout << "Analyze scaffolds" << std::endl;
    for (auto* bam : bamFiles) {
        if (bam && bam->bamType() == "jumps") {
            analyzeStrays(*bam);
        }
    }
}

// =============================================================================
// Circle detection (matching Scala findHgapCircles)
// =============================================================================

std::vector<EndAlignment> Scaffold::findEndAlignments(BamFile& bam) {
    auto refInfo = bam.getRefNamesAndLengths();
    (void)bam.insertSizeMean();     // referenced in Scala but used for logging only
    (void)bam.insertSizeSigma();
    int genomeSize = 0;
    for (const auto& ri : refInfo) genomeSize += ri.second;

    std::cout << "analyzing unpaired in " << bam.path() << std::endl;

    std::vector<EndAlignment> eaList;

    // Iterate reads at contig ends
    for (const auto& ri : refInfo) {
        const std::string& cName = ri.first;
        int cLen = ri.second;
        int endLen = std::min(EndAlignment::EndMaxLength, cLen / 2);

        // Left end: [0, endLen)
        Region leftEnd(cName, 0, endLen);
        auto leftReads = bam.readsInRegion(leftEnd);
        for (const auto& read : leftReads) {
            if (read.unmapped) continue;
            EndAlignment ea(read, cName, cLen);
            if (ea.cantilevered()) eaList.push_back(ea);
        }

        // Right end: [cLen - endLen, cLen)
        Region rightEnd(cName, std::max(0, cLen - endLen), cLen);
        auto rightReads = bam.readsInRegion(rightEnd);
        for (const auto& read : rightReads) {
            if (read.unmapped) continue;
            EndAlignment ea(read, cName, cLen);
            if (ea.cantilevered()) eaList.push_back(ea);
        }
    }

    return eaList;
}

std::unordered_map<std::string, int> Scaffold::findCircles(const std::vector<EndAlignment>& eas) {
    // Group by read name
    std::unordered_map<std::string, std::vector<EndAlignment>> byRead;
    for (const auto& ea : eas) {
        byRead[ea.readName].push_back(ea);
    }

    // Group by pair ends
    // Key: (contigEndName1, contigEndName2, rc xor) — use string key
    std::map<std::string, std::vector<EndAlignmentPair>> byPair;

    for (const auto& [name, caList] : byRead) {
        if (caList.size() <= 1) continue;
        for (size_t i = 0; i < caList.size(); i++) {
            for (size_t j = i + 1; j < caList.size(); j++) {
                EndAlignmentPair caPair(caList[i], caList[j]);
                std::string c1 = caPair.ca1.contigEndName();
                std::string c2 = caPair.ca2.contigEndName();
                bool xorRc = caPair.ca1.rc != caPair.ca2.rc;
                std::string key = c1 + "\t" + c2 + "\t" + (xorRc ? "1" : "0");
                byPair[key].push_back(caPair);
            }
        }
    }

    // Sort by pair list length descending
    std::vector<std::pair<std::string, std::vector<EndAlignmentPair>>> byPairSorted;
    for (auto& [key, pairList] : byPair) {
        byPairSorted.emplace_back(key, std::move(pairList));
    }
    std::sort(byPairSorted.begin(), byPairSorted.end(),
              [](const auto& a, const auto& b) { return a.second.size() > b.second.size(); });

    std::vector<std::vector<EndAlignmentPair>> circleLists;
    for (const auto& [pairEnds, pairList] : byPairSorted) {
        if (static_cast<int>(pairList.size()) < EndMinLinks) continue;
        const auto& pair = pairList.front();
        if (pair.ca1.contigName == pair.ca2.contigName &&
            pair.ca1.contigEndName() != pair.ca2.contigEndName()) {
            circleLists.push_back(pairList);
        }
    }

    std::unordered_map<std::string, int> circles;
    for (const auto& circleList : circleLists) {
        std::string contig = circleList.front().ca1.contigName;
        std::vector<EndAlignmentPair> sortedList = circleList;
        std::sort(sortedList.begin(), sortedList.end(),
                  [](const EndAlignmentPair& a, const EndAlignmentPair& b) {
                      return a.impliedLength() < b.impliedLength();
                  });

        std::vector<int> lengths;
        for (const auto& pair : sortedList) lengths.push_back(pair.impliedLength());
        int estimatedLength = estimateLength(lengths);
        circles[contig] = estimatedLength;

        if (Pilon::verbose) {
            std::cout << sortedList[0].ca1.contigEndName() << " "
                      << sortedList[0].ca2.contigEndName() << " "
                      << sortedList.size() << " " << estimatedLength << std::endl;
            for (const auto& pair : sortedList) {
                std::cout << "  " << pair.impliedLength() << std::endl;
            }
        }
    }

    if (!circles.empty()) {
        std::cout << "Candidate circles: ";
        for (const auto& [name, len] : circles)
            std::cout << name << ":" << len << " ";
        std::cout << std::endl;
    }

    return circles;
}

std::unordered_map<std::string, int> Scaffold::findHgapCircles(std::vector<BamFile*>& bamFiles) {
    if (!Pilon::fixCircles) return {};
    std::vector<EndAlignment> allEas;
    for (auto* bam : bamFiles) {
        if (bam && bam->bamType() == "unpaired") {
            auto eas = findEndAlignments(*bam);
            allEas.insert(allEas.end(), eas.begin(), eas.end());
        }
    }
    return findCircles(allEas);
}

int Scaffold::estimateLength(const std::vector<int>& estimates) {
    int start = 0;
    int n = static_cast<int>(estimates.size());
    while (start < n) {
        int estimate = estimates[start];
        for (int i = start + 1; i < n; i++) {
            if (estimates[i] > estimate + estimate / 5 || i == n - 1) {
                if (i + 1 - start >= EndMinLinks)
                    return estimates[(start + i) / 2];
                else {
                    start = i;
                    break;
                }
            }
        }
    }
    return estimates[n / 2];
}

} // namespace pilon
