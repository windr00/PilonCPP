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

#include "assembler.h"
#include "bases.h"
#include "pilon.h"
#include <algorithm>
#include <iostream>
#include <sstream>

namespace pilon {

int Assembler::K = 47;

Assembler::Assembler(int minDepth_)
    : minDepth_(minDepth_), nReads_(0), nBases_(0), loopLength_(0) {}

void Assembler::addReads(const std::vector<BamRead>& reads) {
    for (const auto& read : reads) {
        addRead(read);
    }
}

void Assembler::addPair(const BamRead& r1, const BamRead& r2) {
    addRead(r1);
    addRead(r2);
}

void Assembler::addRead(const BamRead& r) {
    const std::string& bases = r.bases;
    int length = static_cast<int>(bases.size());
    
    if (length > K) {
        std::vector<uint8_t> quals = r.quals;
        if (quals.empty()) {
            quals.assign(length, Pilon::defaultQual);
        }
        
        int mq = r.mappingQuality;
        addToPileups(bases, quals, mq);

        std::string rcBases = Bases::reverseComplement(bases);
        std::vector<uint8_t> rcQuals = quals;
        std::reverse(rcQuals.begin(), rcQuals.end());
        addToPileups(rcBases, rcQuals, mq);

        nReads_++;
        nBases_ += bases.length();
        
        if (Pilon::verbose && nReads_ % 10000 == 0) {
            std::cout << "..." << nReads_ << std::flush;
        }
        if (nReads_ % 100000 == 0) {
            prunePileups(2);
        }
    }
}

void Assembler::addToPileups(const std::string& bases, const std::vector<uint8_t>& quals, int mq) {
    int length = static_cast<int>(bases.length());
    for (int offset = 0; offset <= length - K - 1; offset++) {
        std::string kmer = bases.substr(offset, K);
        if (pileups.find(kmer) == pileups.end()) {
            pileups[kmer] = PileUp();
        }
        pileups[kmer].add(bases[offset + K], quals[offset + K], mq);
    }
}

void Assembler::addSeq(const std::string& bases) {
    int mq = 10;
    std::vector<uint8_t> quals(bases.length(), 10);
    addToPileups(bases, quals, mq);
    std::string rcBases = Bases::reverseComplement(bases);
    addToPileups(rcBases, quals, mq);
}

void Assembler::addGraphSeq(const std::string& bases) {
    graphSeq(bases);
    std::string rcBases = Bases::reverseComplement(bases);
    graphSeq(rcBases);
}

void Assembler::addGraphSeqs(const std::vector<std::string>& seqs) {
    for (const auto& seq : seqs) {
        addGraphSeq(seq);
    }
}

void Assembler::graphSeq(const std::string& bases) {
    int length = static_cast<int>(bases.length());
    for (int offset = 0; offset <= length - K - 1; offset++) {
        std::string k = bases.substr(offset, K);
        std::string nextK = k.substr(1) + bases[offset + K];

        auto it = kGraph.find(k);
        if (it != kGraph.end() && it->second != nextK) {
            addLink(altGraph, k, nextK, 1);
        } else {
            addLink(kGraph, k, nextK, 1);
        }
    }
}

void Assembler::addLink(KmerGraph& g, const Kmer& k1, const Kmer& k2, int) {
    g[k1] = k2;
}

void Assembler::buildGraph() {
    if (Pilon::debug) {
        std::cout << "building kmer Graph" << std::endl;
    }

    for (auto& kv : pileups) {
        const std::string& k = kv.first;
        PileUp& pu = kv.second;
        
        if (pu.depth() >= minDepth_) {
            auto bc = pu.baseCall();
            std::string prefix = k.substr(1);
            std::string nextK = prefix + bc.base;
            int weight = static_cast<int>(pu.baseCount.sums[bc.baseIndex]);
            
            auto it = kGraph.find(k);
            if (it != kGraph.end() && it->second != nextK) {
                addLink(altGraph, k, nextK, weight);
            } else {
                addLink(kGraph, k, nextK, weight);
            }
            
            if (!bc.homo) {
                addLink(altGraph, k, prefix + bc.altBase, 
                       static_cast<int>(pu.baseCount.sums[bc.altBaseIndex]));
            }
        }
    }
    
    // Free memory
    pileups.clear();

    if (Pilon::debug) {
        std::cout << "kmer graph: t=" << kGraph.size() 
                  << " mt=" << altGraph.size() << std::endl;
    }
}

void Assembler::prunePileups(int minCount) {
    if (Pilon::debug) {
        std::cout << "[prune " << pileups.size();
    }
    
    auto it = pileups.begin();
    while (it != pileups.end()) {
        if (it->second.count() < minCount) {
            it = pileups.erase(it);
        } else {
            ++it;
        }
    }
    
    if (Pilon::debug) {
        std::cout << "->" << pileups.size() << "]" << std::flush;
    }
}

std::string Assembler::kmerPathString(const std::vector<std::string>& kmers, bool prependLength) {
    if (kmers.empty()) return "";
    
    std::string pathStr = kmers.front();
    for (size_t i = 1; i < kmers.size(); i++) {
        pathStr += kmers[i].back();
    }
    
    if (prependLength) {
        return "(" + std::to_string(pathStr.length()) + ")" + pathStr;
    }
    return pathStr;
}

void Assembler::noteKmerLoop(int loopIndex, const std::vector<std::string>& kmers) {
    int length = loopIndex + 1;
    if (loopLength_ == 0 || length < loopLength_) {
        loopLength_ = length;
        std::vector<std::string> subset(kmers.begin(), kmers.begin() + loopIndex + 1);
        loopSequence_ = kmerPathString(subset).substr(0, length);
        
        if (Pilon::verbose) {
            std::cout << "# loop " << loopLength_ << ": " << loopSequence_ << std::endl;
        }
    }
}

std::vector<std::vector<std::string>>
Assembler::kmerPathsForward(std::vector<std::string> kmersIn, int branches) {
    while (true) {
        const std::string& kmer = kmersIn.front();
        
        auto it = kGraph.find(kmer);
        if (it == kGraph.end()) {
            return {kmersIn};
        }
        
        auto altIt = altGraph.find(kmer);
        if (altIt != altGraph.end()) {
            std::string next1 = it->second;
            std::string next2 = altIt->second;
            
            bool seen1 = std::find(kmersIn.begin() + 1, kmersIn.end(), next1) != kmersIn.end();
            bool seen2 = std::find(kmersIn.begin() + 1, kmersIn.end(), next2) != kmersIn.end();

            if (seen1 || seen2) {
                int loop = 0;
                auto pos1 = std::find(kmersIn.begin(), kmersIn.end(), next1);
                auto pos2 = std::find(kmersIn.begin(), kmersIn.end(), next2);
                if (pos1 != kmersIn.end()) loop = std::max(loop, static_cast<int>(pos1 - kmersIn.begin()));
                if (pos2 != kmersIn.end()) loop = std::max(loop, static_cast<int>(pos2 - kmersIn.begin()));
                
                noteKmerLoop(loop, kmersIn);
                if (Pilon::trSafe) return {kmersIn};
            }

            if (seen1 && seen2) return {kmersIn};
            else if (seen1 && !seen2) kmersIn.insert(kmersIn.begin(), next2);
            else if (seen2 && !seen1) kmersIn.insert(kmersIn.begin(), next1);
            else {
                if (branches < maxBranches) {
                    std::vector<std::string> path1 = {next1};
                    path1.insert(path1.end(), kmersIn.begin(), kmersIn.end());
                    auto paths1 = kmerPathsForward(path1, branches + 1);
                    std::vector<std::string> path2 = {next2};
                    path2.insert(path2.end(), kmersIn.begin(), kmersIn.end());
                    auto paths2 = kmerPathsForward(path2, branches + 1);
                    paths1.insert(paths1.end(), paths2.begin(), paths2.end());
                    return paths1;
                }
                return {kmersIn};
            }
        } else {
            std::string next = it->second;
            auto pos = std::find(kmersIn.begin(), kmersIn.end(), next);
            if (pos != kmersIn.end()) {
                noteKmerLoop(static_cast<int>(pos - kmersIn.begin()), kmersIn);
                if (Pilon::trSafe) return {kmersIn};

                int nextCount = std::count(kmersIn.begin(), kmersIn.end(), next);
                if (nextCount > 1) return {kmersIn};
            }
            kmersIn.insert(kmersIn.begin(), next);
        }
    }
    return {};
}

std::vector<std::string> Assembler::pathsForward(const std::string& startingKmer) {
    if (startingKmer.length() != static_cast<size_t>(K)) {
        std::cerr << "Error: starting kmer must be size K" << std::endl;
        return {};
    }
    
    if (kGraph.empty()) buildGraph();
    
    if (Pilon::debug) {
        std::cout << "pathsForward: " << startingKmer << std::flush;
    }
    
    auto kmerPaths = kmerPathsForward({startingKmer});
    std::vector<std::string> paths;
    for (const auto& kp : kmerPaths) {
        paths.push_back(kmerPathString(kp));
    }
    std::sort(paths.begin(), paths.end(), 
              [](const std::string& a, const std::string& b) { return a.length() > b.length(); });
    
    if (Pilon::debug) {
        std::cout << ": " << paths.size() << " paths" << std::endl;
        for (const auto& p : paths) {
            std::cout << "  [" << p.length() << "]" << p << std::endl;
        }
    }
    
    return paths;
}

std::vector<std::string> Assembler::pathsReverse(const std::string& startingKmer) {
    auto paths = pathsForward(Bases::reverseComplement(startingKmer));
    for (auto& p : paths) {
        p = Bases::reverseComplement(p);
    }
    return paths;
}

std::vector<std::string> Assembler::tryForward(const std::string& anchor) {
    if (Pilon::debug) {
        std::cout << "tryForward: [" << anchor.length() << "]" << anchor << std::endl;
    }
    
    if (anchor.length() < static_cast<size_t>(K)) {
        return {anchor};
    } else {
        std::string startingKmer = anchor.substr(0, K);
        auto paths = pathsForward(startingKmer);
        
        if (!paths.empty() && paths[0].length() > anchor.length() + minExtend) {
            return paths;
        } else {
            auto subPaths = tryForward(anchor.substr(1));
            for (auto& p : subPaths) {
                p = startingKmer + p;
            }
            return subPaths;
        }
    }
}

std::vector<std::string> Assembler::tryReverse(const std::string& anchor) {
    if (Pilon::debug) {
        std::cout << "tryReverse" << std::endl;
    }
    
    std::string rcAnchor = Bases::reverseComplement(anchor);
    auto paths = tryForward(rcAnchor);
    for (auto& p : paths) {
        p = Bases::reverseComplement(p);
    }
    
    if (Pilon::debug) {
        std::cout << "RC " << paths.size() << " paths" << std::endl;
        for (const auto& p : paths) {
            std::cout << "  [" << p.length() << "]" << p << std::endl;
        }
    }
    
    return paths;
}

std::tuple<std::vector<std::string>, std::vector<std::string>, std::string>
Assembler::multiBridge(const std::string& left, const std::string& right) {
    auto pathsForward = tryForward(left);
    auto pathsReverse = tryReverse(right);
    return {pathsForward, pathsReverse, loopSequence_};
}

std::vector<std::string> Assembler::novel(Assembler& ref) {
    auto novelKmers = [&](const std::string& seq) -> int {
        int count = 0;
        for (size_t i = 0; i + K <= seq.length(); i++) {
            std::string kmer = seq.substr(i, K);
            if (ref.kGraph.find(kmer) == ref.kGraph.end()) {
                count++;
            }
        }
        return count;
    };

    if (Pilon::verbose) {
        std::cout << "Assembling novel sequence" << std::endl;
    }

    prunePileups();
    // Save k-mer list before building graph (buildGraph clears pileups)
    std::vector<std::string> kmers;
    for (const auto& kv : pileups) kmers.push_back(kv.first);
    if (kGraph.empty()) buildGraph();

    std::unordered_set<std::string> usedKmers;
    std::vector<std::string> paths;
    int n = 0;

    for (const auto& kmer : kmers) {
        if (usedKmers.find(kmer) != usedKmers.end()) continue;
        
        auto forwards = pathsForward(kmer);
        auto reverses = pathsReverse(kmer);
        
        std::string forward = forwards.empty() ? "" : forwards.front();
        std::string reverse = reverses.empty() ? "" : reverses.front();
        std::string path = reverse + forward.substr(K);
        
        for (size_t i = 0; i + K <= path.length(); i++) {
            std::string k = path.substr(i, K);
            usedKmers.insert(k);
            usedKmers.insert(Bases::reverseComplement(k));
        }
        
        if (path.length() >= minNovel) {
            paths.push_back(path);
        }
        
        n++;
    }
    
    std::sort(paths.begin(), paths.end(),
              [](const std::string& a, const std::string& b) { return a.length() > b.length(); });
    
    std::vector<std::string> result;
    for (const auto& path : paths) {
        int kLength = static_cast<int>(path.length()) - (K - 1);
        int novel = novelKmers(path);
        int novelPct = kLength > 0 ? (novel * 100) / kLength : 0;
        
        if (novel >= minNovel && novelPct >= minNovelPct) {
            ref.addGraphSeq(path);
            if (Pilon::verbose) {
                std::cout << "novel " << path.length() << " " << novelPct << "% " << path << std::endl;
            }
            result.push_back(path);
        }
    }
    
    return result;
}

std::string Assembler::toString() const {
    std::ostringstream oss;
    oss << "<assembler K=" << K << " nReads=" << nReads_ 
        << " nBases=" << nBases_ << " nKmers=" << pileups.size() << ">";
    return oss.str();
}

} // namespace pilon
