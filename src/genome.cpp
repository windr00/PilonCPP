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

#include "genome.h"
#include "utils.h"
#include "pilon.h"
#include "bamfile.h"
#include "gapfiller.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>
#include <queue>
#include <functional>
#include <vector>
#include <algorithm>
#include <cmath>

namespace pilon {

// Thread-safe output mutex
static std::mutex coutMutex;

// =============================================================================
// Helper: process a single chunk (called from worker threads)
// =============================================================================
static void processChunk(const std::string& name,
                         const std::string& seq,
                         int chunkStart, int chunkStop,
                         std::vector<BamFile*>& bamFiles,
                         std::vector<GenomeRegion>& results,
                         std::mutex& resultsMutex,
                         std::atomic<int>& completedChunks,
                         int totalChunks) {
    std::vector<BamFile*> threadBams;
    for (auto* bam : bamFiles) {
        if (bam) {
            auto* threadBam = new BamFile(bam->path(), bam->bamType(), bam->subType());
            threadBam->open();
            threadBams.push_back(threadBam);
        }
    }

    GenomeRegion region(name, chunkStart, chunkStop,
                        seq.substr(chunkStart, chunkStop - chunkStart),
                        Pilon::minDepth);

    for (auto* bam : threadBams) {
        bam->process(region);
    }

    region.postProcess();

    if (Pilon::fixSnps || Pilon::fixIndels || Pilon::fixGaps || Pilon::fixLocal) {
        region.identifyAndFixIssues();
    }

    {
        std::lock_guard<std::mutex> lock(resultsMutex);
        results.push_back(std::move(region));
    }

    int done = ++completedChunks;
    if (Pilon::verbose || done % 10 == 0 || done == totalChunks) {
        std::lock_guard<std::mutex> lock(coutMutex);
        std::cout << "  [" << done << "/" << totalChunks << "] Chunk "
                  << chunkStart << "-" << chunkStop << " done" << std::endl;
    }

    for (auto* bam : threadBams) {
        bam->close();
        delete bam;
    }
}

// =============================================================================
// GenomeRegion constructor
// =============================================================================
GenomeRegion::GenomeRegion(const std::string& name, int start, int stop,
                           const std::string& bases, double minDepth)
    : name(name), start(start), stop(stop), contigBases(bases),
      originalBases(bases), bases(bases), minDepth(minDepth),
      physCovStart(0), insertSizeStart(0), readCount(0), baseCount(0) {
    int sz = stop - start;
    pileUps.resize(sz);
    confirmed.resize(sz, false);
    ambiguous.resize(sz, false);
    changed.resize(sz, false);
    deleted.resize(sz, false);
    excluded.resize(sz, false);
    coverage_arr.resize(sz, 0);
    badCoverage_arr.resize(sz, 0);
    clips_arr.resize(sz, 0);
    insertSize_arr.resize(sz, 0);
    physCoverage_arr.resize(sz, 0);
    fragCoverage_arr.resize(sz, 0);
    weightedQual_arr.resize(sz, 0);
    weightedMq_arr.resize(sz, 0);
}

// =============================================================================
// Basic accessors
// =============================================================================
int GenomeRegion::size() const { return stop - start; }

char GenomeRegion::baseAt(int pos) const {
    if (pos >= 0 && pos < static_cast<int>(contigBases.size()))
        return contigBases[pos];
    return 'N';
}

char GenomeRegion::originalBaseAt(int pos) const {
    if (pos >= 0 && pos < static_cast<int>(originalBases.size()))
        return originalBases[pos];
    return 'N';
}

char GenomeRegion::refBase(int pos) const {
    int localPos = pos - start;
    return baseAt(localPos);
}

std::string GenomeRegion::subString(int start, int length) const {
    if (start < 0 || start >= static_cast<int>(contigBases.size())) return "";
    int len = std::min(length, static_cast<int>(contigBases.size()) - start);
    return contigBases.substr(start, len);
}

std::string GenomeRegion::refSubString(int start, int length) const {
    return subString(start, length);
}

PileUp& GenomeRegion::pileUpRegion(int index) {
    return pileUps[index];
}

const PileUp& GenomeRegion::pileUpRegion(int index) const {
    return pileUps[index];
}

int GenomeRegion::locus(int index) const {
    return start + index;
}

bool GenomeRegion::inRegion(int locus) const {
    return locus >= start && locus < stop;
}

bool GenomeRegion::beforeRegion(int locus) const {
    return locus < start;
}

bool GenomeRegion::afterRegion(int locus) const {
    return locus >= stop;
}

int GenomeRegion::index(int locus) const {
    return locus - start;
}

// =============================================================================
// Homo run and nanopore exclude (matching Scala)
// =============================================================================
int GenomeRegion::homoRun(int loc) const {
    if (loc < 0 || loc >= static_cast<int>(contigBases.size())) return 0;
    char baseAtLoc = contigBases[loc];
    for (int i = loc + 1; i < static_cast<int>(contigBases.size()); i++) {
        if (contigBases[i] != baseAtLoc) return i - loc;
    }
    return static_cast<int>(contigBases.size()) - loc;
}

bool GenomeRegion::nanoporeExclude(int idx) const {
    return (idx - 2 >= 0 && idx + 2 < static_cast<int>(contigBases.size()) &&
            contigBases[idx - 2] == 'C' &&
            contigBases[idx - 1] == 'C' &&
            contigBases[idx + 1] == 'G' &&
            contigBases[idx + 2] == 'G');
}

void GenomeRegion::excludeMotifs() {
    bool pb = Pilon::pacbio;
    bool nano = Pilon::nanopore;
    bool lr = pb || nano;

    if (!lr) return;

    for (int i = 0; i < static_cast<int>(contigBases.size()); i++) {
        excluded[i] = (homoRun(i) >= 4) || (nano && nanoporeExclude(i));
    }
}

// =============================================================================
// Physical coverage computation (matching Scala PileUpRegion.computePhysCov)
// =============================================================================
void GenomeRegion::computePhysCov() {
    if (pileUps.empty()) return;
    
    pileUps[0].physCov += physCovStart;
    pileUps[0].insertSize += insertSizeStart;
    
    for (int i = 1; i < static_cast<int>(pileUps.size()); i++) {
        pileUps[i].physCov += pileUps[i - 1].physCov;
        pileUps[i].insertSize += pileUps[i - 1].insertSize;
    }
    
    for (int i = 0; i < static_cast<int>(pileUps.size()); i++) {
        if (pileUps[i].physCov > 0) {
            pileUps[i].insertSize /= pileUps[i].physCov;
        }
    }
}

// =============================================================================
// Post-process: collect changes from pileups (matching Scala postProcess)
// =============================================================================
void GenomeRegion::postProcess() {
    computePhysCov();
    
    // Compute meanCoverage
    int meanCoverage = 0;
    if (!pileUps.empty()) {
        long long totalDepth = 0;
        for (const auto& pu : pileUps) {
            totalDepth += pu.depth();
        }
        meanCoverage = static_cast<int>(totalDepth / pileUps.size());
    }
    
    // Compute minDepth (matching Scala postProcess logic)
    minDepth = Pilon::minDepth >= 1 
        ? static_cast<int>(Pilon::minDepth)
        : std::max(static_cast<int>(std::round(Pilon::minDepth * meanCoverage)), Pilon::minMinDepth);
    
    if (readCount == 0) {
        if (Pilon::verbose) {
            std::cerr << "WARNING: No reads processed in region " << name 
                      << ":" << start << "-" << stop << std::endl;
        }
        return;
    }
    
    // excluded motifs for long reads (matching Scala)
    if (Pilon::longread) {
        excludeMotifs();
    }
    
    // Pass 1: pull out values from pileups & call base changes
    bool fixamb = Pilon::iupac || Pilon::fixAmb;
    
    for (int i = 0; i < size(); i++) {
        const PileUp& pu = pileUps[i];
        long long n = pu.depth();
        auto bc = pu.baseCall();
        char b = bc.base;
        bool homo = bc.homo;
        char r = contigBases[i];
        
        // Store summary stats (matching Scala)
        coverage_arr[i] = n;
        badCoverage_arr[i] = pu.badPair;
        physCoverage_arr[i] = pu.physCov;
        insertSize_arr[i] = pu.insertSize;
        weightedQual_arr[i] = static_cast<int8_t>(pu.weightedQual());
        weightedMq_arr[i] = static_cast<int8_t>(pu.weightedMq());
        clips_arr[i] = pu.clips;
        
        if (n >= minDepth && r != 'N' && !deleted[i] && bc.called()) {
            if (homo && b == r && bc.highConfidence() && !bc.indel) {
                confirmed[i] = true;
            } else if (bc.isInsertion() && bc.homoIndel) {
                changes_.push_back({i, INS, pu});
            } else if (bc.isDeletion() && bc.homoIndel) {
                changes_.push_back({i, DEL, pu});
                // Scala: mark deleted positions for downstream positions
                for (int j = 1; j < static_cast<int>(bc.deletion.length()); j++) {
                    int delIdx = i + j;
                    if (delIdx < static_cast<int>(deleted.size())) {
                        deleted[delIdx] = true;
                        pileUps[delIdx].deletions += pu.deletions;
                    }
                }
            } else if (b != r && bc.score > 0) {
                if (homo) {
                    changes_.push_back({i, SNP, pu});
                } else if (fixamb || bc.altBase != r) {
                    changes_.push_back({i, AMB, pu});
                }
            }
        }
    }
    
    // Pass 2: computed values (simplified for now)
    // Scala also computes copyNumber, fragCoverageDist here
}

// =============================================================================
// fixFixList: Sort and remove overlapping fixes (matching Scala)
// =============================================================================
static std::vector<GenomeRegion::Fix> fixFixList(const std::vector<GenomeRegion::Fix>& inList) {
    std::vector<GenomeRegion::Fix> fixes = inList;
    std::sort(fixes.begin(), fixes.end(),
              [](const GenomeRegion::Fix& a, const GenomeRegion::Fix& b) {
                  return std::get<0>(a) < std::get<0>(b);
              });
    
    std::vector<GenomeRegion::Fix> outList;
    
    while (!fixes.empty()) {
        if (fixes.size() >= 2) {
            const auto& fix1 = fixes[0];
            const auto& fix2 = fixes[1];
            
            int region1End = std::get<0>(fix1) + std::max(static_cast<int>(std::get<1>(fix1).length()) - 1, 0);
            int region2End = std::get<0>(fix2) + std::max(static_cast<int>(std::get<1>(fix2).length()) - 1, 0);
            
            if (std::get<0>(fix1) <= region2End && std::get<0>(fix2) <= region1End) {
                int fix1len = std::get<1>(fix1).length() + std::get<2>(fix1).length();
                int fix2len = std::get<1>(fix2).length() + std::get<2>(fix2).length();
                
                if (fix1len >= fix2len) {
                    // Scala: fix1 wins, drop fix2, keep fix1 for further comparison
                    fixes.erase(fixes.begin() + 1);
                } else {
                    outList.push_back(fix2);
                    fixes.erase(fixes.begin(), fixes.begin() + 2);
                }
            } else {
                outList.push_back(fix1);
                fixes.erase(fixes.begin());
            }
        } else {
            outList.push_back(fixes[0]);
            fixes.clear();
        }
    }
    
    return outList;
}

// =============================================================================
// fixIssues: Apply fixes to bases (matching Scala fixIssues)
// =============================================================================
void GenomeRegion::identifyAndFixIssues() {
    std::vector<Fix> snpFixList;
    std::vector<Fix> smallFixList;
    std::vector<Fix> bigFixList;
    
    if (Pilon::verbose) {
        std::cout << "  Changes found: " << changes_.size() << std::endl;
    }
    
    for (const auto& change : changes_) {
        int i = change.index;
        ChangeKind kind = change.kind;
        const PileUp& pu = change.pu;
        
        int loc = locus(i);
        char rBase = contigBases[i];
        auto bc = pu.baseCall();
        char cBase = bc.base;
        
        if (!excluded[i]) {
            switch (kind) {
                case SNP:
                    if (Pilon::fixSnps) {
                        snpFixList.push_back({loc, std::string(1, rBase), std::string(1, cBase)});
                    }
                    break;
                case AMB:
                    if (Pilon::fixSnps && !Pilon::longread) {
                        if (Pilon::iupac) {
                            smallFixList.push_back({loc, std::string(1, rBase),
                                                    std::string(1, bc.iupacBase())});
                        } else {
                            snpFixList.push_back({loc, std::string(1, rBase), std::string(1, cBase)});
                        }
                    }
                    break;
                case INS:
                    if (Pilon::fixIndels) {
                        smallFixList.push_back({loc, "", bc.insertion});
                    }
                    break;
                case DEL:
                    if (Pilon::fixIndels) {
                        smallFixList.push_back({loc, bc.deletion, ""});
                    }
                    break;
            }
        }
    }
    
    // Report stats (matching Scala)
    int ins = 0, dels = 0;
    int insBases = 0, delBases = 0;
    for (const auto& fix : smallFixList) {
        if (std::get<1>(fix).empty() && !std::get<2>(fix).empty()) {
            ins++;
            insBases += std::get<2>(fix).length();
        } else if (!std::get<1>(fix).empty() && std::get<2>(fix).empty()) {
            dels++;
            delBases += std::get<1>(fix).length();
        }
    }
    
    // Apply SNP fixes first, then small+big fixes (matching Scala)
    if (!snpFixList.empty()) {
        auto snpFixes = fixFixList(snpFixList);
        for (auto it = snpFixes.rbegin(); it != snpFixes.rend(); ++it) {
            int locus_global = std::get<0>(*it);
            const std::string& was = std::get<1>(*it);
            const std::string& patch = std::get<2>(*it);
            int startIdx = locus_global - start;
            for (size_t i = 0; i < was.length(); i++) {
                int iNew = startIdx + i;
                if (iNew >= 0 && iNew < static_cast<int>(bases.length())) {
                    bases[iNew] = patch[i];
                }
            }
        }
    }
    
    if (!smallFixList.empty() || !bigFixList.empty()) {
        std::vector<Fix> allFixes;
        allFixes.insert(allFixes.end(), smallFixList.begin(), smallFixList.end());
        allFixes.insert(allFixes.end(), bigFixList.begin(), bigFixList.end());
        
        auto sortedFixes = fixFixList(allFixes);
        for (auto it = sortedFixes.rbegin(); it != sortedFixes.rend(); ++it) {
            int locus_global = std::get<0>(*it);
            const std::string& was = std::get<1>(*it);
            const std::string& patch = std::get<2>(*it);
            int startIdx = locus_global - start;
            
            if (was.length() == patch.length()) {
                for (size_t i = 0; i < was.length(); i++) {
                    int iNew = startIdx + i;
                    if (iNew >= 0 && iNew < static_cast<int>(bases.length())) {
                        bases[iNew] = patch[i];
                    }
                }
            } else {
                if (startIdx >= 0 && startIdx <= static_cast<int>(bases.length())) {
                    std::string before = bases.substr(0, startIdx);
                    std::string after = (startIdx + static_cast<int>(was.length()) <= static_cast<int>(bases.length()))
                        ? bases.substr(startIdx + was.length())
                        : "";
                    bases = before + patch + after;
                }
            }
        }
    }
    
    // Store fixes for VCF output
    fixes = snpFixList;
    fixes.insert(fixes.end(), smallFixList.begin(), smallFixList.end());
    fixes.insert(fixes.end(), bigFixList.begin(), bigFixList.end());
}

// =============================================================================
// Write VCF record
// =============================================================================
void GenomeRegion::writeVcf(FILE* writer) const {
    if (!writer) return;
    for (int i = 0; i < size(); i++) {
        const PileUp& pu = pileUps[i];
        auto bc = pu.baseCall();
        if (!bc.called()) continue;
        if (pu.depth() < static_cast<long long>(minDepth)) continue;
        
        char refB = contigBases[i];
        int pos = locus(i) + 1; // VCF is 1-based
        
        if (Pilon::fixSnps && bc.called() && !bc.isInsertion() && !bc.isDeletion()) {
            if (!bc.baseMatch(refB)) {
                std::string refStr(1, refB);
                std::string altStr(1, bc.base);
                fprintf(writer, "%s\t%d\t.\t%s\t%s\t.\tPASS\tDP=%lld\n",
                        name.c_str(), pos, refStr.c_str(), altStr.c_str(),
                        pu.depth());
            }
        }
        
        if (Pilon::fixIndels) {
            if (bc.isInsertion()) {
                std::string refStr(1, refB);
                std::string altStr = refB + bc.insertion;
                fprintf(writer, "%s\t%d\t.\t%s\t%s\t.\tPASS\tDP=%lld\n",
                        name.c_str(), pos, refStr.c_str(), altStr.c_str(),
                        pu.depth());
            } else if (bc.isDeletion()) {
                std::string refStr = refB + bc.deletion;
                std::string altStr(1, refB);
                fprintf(writer, "%s\t%d\t.\t%s\t%s\t.\tPASS\tDP=%lld\n",
                        name.c_str(), pos, refStr.c_str(), altStr.c_str(),
                        pu.depth());
            }
        }
    }
}

// =============================================================================
// Write changes
// =============================================================================
void GenomeRegion::writeChanges(FILE* writer, const std::string& newName, int& offset) const {
    for (const auto& fix : fixes) {
        int pos = std::get<0>(fix) + offset + 1; // 1-based
        const std::string& refSeq = std::get<1>(fix);
        const std::string& altSeq = std::get<2>(fix);
        fprintf(writer, "%s\t%d\t%s\t%s\n",
                newName.c_str(), pos, refSeq.c_str(), altSeq.c_str());
    }
    offset += static_cast<int>(bases.size()) - size();
}

// =============================================================================
// GenomeFile implementation
// =============================================================================
GenomeFile::GenomeFile(const std::string& path, const std::string& targets)
    : targets_(targets) {
    parseFasta(path);
}

void GenomeFile::parseFasta(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open genome file: " << path << std::endl;
        exit(1);
    }

    std::string line;
    std::string currentName;
    std::string currentSeq;

    while (std::getline(file, line)) {
        line = Utils::trim(line);
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!currentName.empty()) {
                contigs_.push_back({currentName, currentSeq});
                currentSeq.clear();
            }
            std::istringstream iss(line.substr(1));
            iss >> currentName;
        } else {
            for (char c : line) {
                currentSeq += std::toupper(c);
            }
        }
    }
    if (!currentName.empty()) {
        contigs_.push_back({currentName, currentSeq});
    }
}

std::vector<std::pair<std::string, std::string>> GenomeFile::getContigs() const {
    return contigs_;
}

std::vector<std::pair<std::string, int>> GenomeFile::getContigSizes() const {
    std::vector<std::pair<std::string, int>> sizes;
    for (const auto& contig : contigs_) {
        sizes.push_back({contig.first, static_cast<int>(contig.second.size())});
    }
    return sizes;
}

std::vector<std::pair<std::string, std::string>> GenomeFile::parseTargets() {
    if (targets_.empty()) return contigs_;
    std::string content = Utils::readFile(targets_);
    std::istringstream stream(content);
    std::vector<std::pair<std::string, std::string>> result;
    std::string line;
    while (std::getline(stream, line)) {
        line = Utils::trim(line);
        if (line.empty() || line[0] == '#') continue;
        auto parts = Utils::split(line, '\t');
        if (parts.size() >= 1) {
            std::string name = Utils::trim(parts[0]);
            for (const auto& contig : contigs_) {
                if (contig.first == name) {
                    result.push_back(contig);
                    break;
                }
            }
        }
    }
    return result;
}

void GenomeFile::processRegions(std::vector<BamFile*>& bamFiles) {
    auto contigs = targets_.empty() ? contigs_ : parseTargets();
    int numThreads = Pilon::threads;

    if (numThreads <= 1) {
        for (const auto& contig : contigs) {
            const std::string& name = contig.first;
            const std::string& seq = contig.second;
            int length = static_cast<int>(seq.size());
            std::cout << "Processing " << name << " (" << length << " bp)" << std::endl;

            for (int chunkStart = 0; chunkStart < length; chunkStart += Pilon::chunkSize) {
                int chunkStop = std::min(chunkStart + Pilon::chunkSize, length);
                GenomeRegion region(name, chunkStart, chunkStop,
                                   seq.substr(chunkStart, chunkStop - chunkStart),
                                   Pilon::minDepth);

                for (auto* bam : bamFiles) {
                    if (bam) bam->process(region);
                }

                region.postProcess();

                if (Pilon::fixSnps || Pilon::fixIndels || Pilon::fixGaps || Pilon::fixLocal) {
                    region.identifyAndFixIssues();
                }

                processedRegions_.push_back(std::move(region));

                if (Pilon::verbose) {
                    std::cout << "  Chunk " << chunkStart << "-" << chunkStop << " done" << std::endl;
                }
            }
        }
        return;
    }

    // Multi-threaded mode
    std::cout << "Using " << numThreads << " threads for processing" << std::endl;
    struct ChunkTask { std::string name; std::string seq; int chunkStart; int chunkStop; };
    std::vector<ChunkTask> tasks;
    for (const auto& contig : contigs) {
        const std::string& name = contig.first;
        const std::string& seq = contig.second;
        int length = static_cast<int>(seq.size());
        std::cout << "Processing " << name << " (" << length << " bp)" << std::endl;
        for (int chunkStart = 0; chunkStart < length; chunkStart += Pilon::chunkSize) {
            int chunkStop = std::min(chunkStart + Pilon::chunkSize, length);
            tasks.push_back({name, seq, chunkStart, chunkStop});
        }
    }

    int totalChunks = static_cast<int>(tasks.size());
    std::cout << "Total chunks: " << totalChunks << std::endl;
    std::mutex resultsMutex;
    std::atomic<int> completedChunks(0);
    processedRegions_.reserve(totalChunks);
    std::vector<std::thread> workers;
    workers.reserve(numThreads);
    std::mutex taskMutex;
    std::queue<size_t> taskQueue;
    for (size_t i = 0; i < tasks.size(); i++) taskQueue.push(i);

    auto workerFunc = [&]() {
        while (true) {
            size_t taskIdx;
            {
                std::lock_guard<std::mutex> lock(taskMutex);
                if (taskQueue.empty()) return;
                taskIdx = taskQueue.front();
                taskQueue.pop();
            }
            const auto& task = tasks[taskIdx];
            processChunk(task.name, task.seq, task.chunkStart, task.chunkStop,
                        bamFiles, processedRegions_, resultsMutex,
                        completedChunks, totalChunks);
        }
    };

    for (int i = 0; i < numThreads; i++) workers.emplace_back(workerFunc);
    for (auto& t : workers) t.join();

    std::sort(processedRegions_.begin(), processedRegions_.end(),
              [](const GenomeRegion& a, const GenomeRegion& b) {
                  if (a.name != b.name) return a.name < b.name;
                  return a.start < b.start;
              });
    std::cout << "Multi-threaded processing complete: " << totalChunks << " chunks" << std::endl;
}

} // namespace pilon
