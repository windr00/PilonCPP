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

int GapFiller::k = 0;

static bool substrEq(const std::string& a, int aOff, const std::string& b, int bOff, int len) {
    if (aOff + len > (int)a.length() || bOff + len > (int)b.length()) return false;
    for (int i = 0; i < len; i++)
        if (a[aOff + i] != b[bOff + i]) return false;
    return true;
}

static std::string findProperOverlap(const std::string& left, const std::string& right, int minOverlap) {
    int ll = (int)left.length();
    int rl = (int)right.length();
    if (ll < minOverlap || rl < minOverlap) return "";
    int maxOverlap = std::min(ll, rl);
    for (int overlap = minOverlap; overlap <= maxOverlap; overlap++)
        if (substrEq(left, ll - overlap, right, 0, overlap))
            return left + right.substr(overlap);
    return "";
}

void GapFiller::fixBreak(GenomeRegion& region, const Region& breakRegion) {
    if (Pilon::verbose)
        std::cout << "  fixBreak at " << breakRegion.toString() << std::endl;

    int minRadius = 3 * Assembler::K;
    int radius = minRadius;

    int start = std::max(breakRegion.start - radius, 0);
    int stop = std::min(breakRegion.stop + radius, region.size());
    if (start >= stop) return;

    std::string leftFlank = region.subString(start, breakRegion.start - start);
    std::string rightFlank = region.subString(breakRegion.stop, stop - breakRegion.stop);
    if (leftFlank.empty() || rightFlank.empty()) return;

    if (Pilon::verbose)
        std::cout << "    radius=" << radius << " left=" << leftFlank.length()
                  << " right=" << rightFlank.length() << std::endl;

    // Recruit reads flanking the break from all BAMs
    Assembler assembler;
    assembler.addSeq(leftFlank);
    assembler.addSeq(rightFlank);

    for (auto* bam : Pilon::bamFiles) {
        if (!bam) continue;
        std::vector<BamRead> flankReads = bam->recruitFlankReads(
            Region(breakRegion.name, start, stop));
        for (auto& br : flankReads) {
            if (!br.bases.empty())
                assembler.addGraphSeq(br.bases.c_str());
        }
        if (Pilon::verbose)
            std::cout << "    recruited " << flankReads.size() << " reads from " << bam->path().c_str() << std::endl;
    }

    k = 2 * Assembler::K + 1;
    int minOverlap = std::min(k, std::min((int)leftFlank.length(), (int)rightFlank.length()));
    if (minOverlap < 10) minOverlap = 10;

    // Try multiBridge using assembled reads
    auto [fwd, rev, loop] = assembler.multiBridge(leftFlank, rightFlank);
    bool bridgeFound = false;
    std::string patch;

    if (!fwd.empty() && !rev.empty()) {
        if (Pilon::verbose)
            std::cout << "    multiBridge: fwd[0]=" << fwd[0].length()
                      << " rev[0]=" << rev[0].length()
                      << " loop=" << loop.length() << std::endl;

        // Try connecting forward path to right flank
        int fwdOverlap = std::min((int)fwd[0].length(), (int)rightFlank.length());
        for (int ov = std::min(fwdOverlap, 200); ov >= minOverlap; ov--) {
            if (substrEq(fwd[0], fwd[0].length() - ov, rightFlank, 0, ov)) {
                patch = fwd[0] + rightFlank.substr(ov);
                bridgeFound = true;
                break;
            }
        }

        // Try connecting reverse path to left flank
        if (!bridgeFound) {
            int revOverlap = std::min((int)rev[0].length(), (int)leftFlank.length());
            for (int ov = std::min(revOverlap, 200); ov >= minOverlap; ov--) {
                if (substrEq(leftFlank, leftFlank.length() - ov, rev[0], 0, ov)) {
                    patch = leftFlank + rev[0].substr(ov);
                    bridgeFound = true;
                    break;
                }
            }
        }

        // Try consensus between forward and reverse paths
        if (!bridgeFound) {
            int minLen = std::min((int)fwd[0].length(), (int)rev[0].length());
            for (int ov = std::min(minLen, 200); ov >= minOverlap; ov--) {
                if (substrEq(fwd[0], fwd[0].length() - ov, rev[0], 0, ov)) {
                    patch = fwd[0] + rev[0].substr(ov);
                    bridgeFound = true;
                    break;
                }
            }
        }

        if (!bridgeFound) {
            std::string combined = leftFlank + fwd[0] + rev[0];
            patch = findProperOverlap(combined, rightFlank, minOverlap);
            if (!patch.empty()) bridgeFound = true;
        }
    }

    if (!bridgeFound && !fwd.empty() && !rev.empty()) {
        patch = findProperOverlap(fwd[0], rev[0], minOverlap);
        if (!patch.empty()) bridgeFound = true;
    }

    if (!bridgeFound && !loop.empty() && loop.length() >= (size_t)minOverlap) {
        if (Pilon::verbose)
            std::cout << "    using loop: " << loop.substr(0, 40) << std::endl;
        patch = loop;
        bridgeFound = true;
    }

    // Fallback: direct properOverlap
    if (!bridgeFound) {
        patch = findProperOverlap(leftFlank, rightFlank, minOverlap);
        if (!patch.empty()) bridgeFound = true;
    }

    if (!bridgeFound) {
        std::string rcRight = rightFlank;
        std::reverse(rcRight.begin(), rcRight.end());
        for (auto& c : rcRight) {
            switch (c) {
                case 'A': c = 'T'; break;
                case 'T': c = 'A'; break;
                case 'C': c = 'G'; break;
                case 'G': c = 'C'; break;
            }
        }
        patch = findProperOverlap(leftFlank, rcRight, minOverlap);
        if (!patch.empty()) bridgeFound = true;
    }

    if (!bridgeFound && !fwd.empty()) {
        patch = findProperOverlap(leftFlank, fwd[0], minOverlap);
        if (!patch.empty()) bridgeFound = true;
    }

    if (bridgeFound) {
        int trimStart = 0;
        while (trimStart < (int)patch.length() && trimStart < radius &&
               start + trimStart < region.stop &&
               region.baseAt(start + trimStart) == patch[trimStart])
            trimStart++;

        int trimEnd = 0;
        int plen = (int)patch.length();
        while (trimEnd < plen && trimEnd < radius &&
               stop - trimEnd - 1 >= 0 &&
               region.baseAt(stop - trimEnd - 1) == patch[plen - 1 - trimEnd])
            trimEnd++;

        std::string trimmedPatch = patch.substr(trimStart, plen - trimStart - trimEnd);
        if ((int)trimmedPatch.length() > 5) {
            if (Pilon::verbose)
                std::cout << "    bridge found: len=" << trimmedPatch.length()
                          << " " << trimmedPatch.substr(0, 40) << "..." << std::endl;
            region.fixBreakRegion(breakRegion, trimmedPatch);
        } else if (Pilon::verbose) {
            std::cout << "    patch too short after trim" << std::endl;
        }
    } else if (Pilon::verbose) {
        std::cout << "    no bridge found" << std::endl;
    }
}

void GapFiller::fixGap(GenomeRegion& region, const Region& gapRegion) {
    fixBreak(region, gapRegion);
}

bool GapFiller::fixCircles(GenomeRegion& region, int gapMargin) {
    if (Pilon::verbose) std::cout << "  fixCircles on " << region.name << std::endl;
    int contigLen = region.size();
    if (contigLen < 1000) return false;

    int endLen = std::min(5000, contigLen / 2);

    struct EndAlign {
        std::string qname;
        int start, end;
        bool rc;
        bool leftCant, rightCant;
    };

    std::vector<EndAlign> endAligns;

    for (auto* bam : Pilon::bamFiles) {
        if (!bam) continue;
        std::vector<BamRead> reads = bam->readsInRegion(Region(region.name, 0, contigLen));
        for (auto& br : reads) {
            if (br.bases.empty()) continue;
            int aStart = br.alignmentStart;
            int aEnd = br.alignmentEnd;

            if (aStart <= endLen && aStart < 100) {
                if (br.negativeStrand)
                    endAligns.push_back({br.readName, aStart, aEnd, br.negativeStrand, true, false});
            }
            if (aEnd >= contigLen - endLen && aEnd > contigLen - 100) {
                if (!br.negativeStrand)
                    endAligns.push_back({br.readName, aStart, aEnd, br.negativeStrand, false, true});
            }
        }
        if (Pilon::verbose)
            std::cout << "    " << reads.size() << " reads, " << endAligns.size()
                      << " cantalevering from " << bam->path().c_str() << std::endl;
    }

    std::map<std::pair<std::string, bool>, EndAlign> pairedEnds;
    bool circularFound = false;
    for (auto& ea : endAligns) {
        auto key = std::make_pair(ea.qname, ea.rc);
        auto it = pairedEnds.find(key);
        if (it != pairedEnds.end()) {
            if ((it->second.leftCant && ea.rightCant) || (it->second.rightCant && ea.leftCant)) {
                if (Pilon::verbose)
                    std::cout << "    circular bridge: " << ea.qname << std::endl;
                circularFound = true;
            }
        } else {
            pairedEnds[key] = ea;
        }
    }

    if (Pilon::verbose) {
        int checkLen = std::min(500, contigLen / 4);
        long leftCov = 0, rightCov = 0;
        int n = 0;
        for (int i = 0; i < checkLen && i < (int)region.coverage_arr.size(); i++) {
            leftCov += region.coverage_arr[i];
            rightCov += region.coverage_arr[region.coverage_arr.size() - 1 - i];
            n++;
        }
        if (n > 0)
            std::cout << "    end_cov avg=" << (leftCov + rightCov) / (2 * n) << std::endl;
        std::cout << "    circular=" << (circularFound ? "yes" : "no") << std::endl;
    }

    return circularFound;
}

void GapFiller::fixNovel(GenomeFile* genome, const std::vector<BamFile*>& bamFiles) {
    if (Pilon::verbose)
        std::cout << "  fixNovel: assembling novel contigs from unmapped reads" << std::endl;

    // Build reference genome graph (matching Scala: genomeGraph = new Assembler(minDepth=1))
    Assembler genomeGraph;
    if (genome) {
        for (auto& contig : genome->getContigs()) {
            genomeGraph.addGraphSeq(contig.second.c_str());
        }
    }

    // Collect unmapped reads from non-jump BAMs (frags and unpaired)
    std::vector<BamRead> unaligned;
    for (auto* bam : Pilon::bamFiles) {
        if (!bam) continue;
        // Scala: skip jump BAMs for novel assembly
        if (bam->bamType() == "jumps") continue;

        std::vector<BamRead> ua = bam->getUnalignedReads();
        unaligned.insert(unaligned.end(), ua.begin(), ua.end());
        if (Pilon::verbose)
            std::cout << "    " << ua.size() << " unaligned from " << bam->path().c_str() << std::endl;
    }

    if (unaligned.empty()) {
        if (Pilon::verbose)
            std::cout << "    no unaligned reads found" << std::endl;
        return;
    }

    // Build assembler from unaligned reads (matching Scala: novelGraph = new Assembler)
    Assembler novelGraph(Pilon::minDepth);  // Use proper depth threshold
    novelGraph.addReads(unaligned);

    // Extract novel contigs by subtracting reference k-mers
    std::vector<std::string> novelContigs = novelGraph.novel(genomeGraph);

    if (Pilon::verbose)
        std::cout << "    found " << novelContigs.size() << " novel contigs" << std::endl;

    // Output novel contigs (matching Scala: write to <output>.novel.fasta)
    std::string novelPath = Pilon::outputFile("novel.fasta");
    FILE* out = fopen(novelPath.c_str(), "w");
    if (out) {
        int count = 0;
        for (auto& nc : novelContigs) {
            if ((int)nc.length() >= 200) {  // Filter: min length 200
                count++;
                fprintf(out, ">novel_contig_%d length=%zu coverage=%.1f\n",
                        count, nc.length(), 1.0);
                for (size_t i = 0; i < nc.length(); i += 80)
                    fprintf(out, "%.*s\n", std::min(80, (int)(nc.length() - i)), nc.c_str() + i);
            }
        }
        fclose(out);
        if (Pilon::verbose)
            std::cout << "    wrote " << count << " novel contigs to " << novelPath << std::endl;
    }
}

} // namespace pilon
