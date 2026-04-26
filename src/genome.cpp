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

namespace pilon {

// Thread-safe output mutex
static std::mutex coutMutex;

// Helper: process a single chunk (called from worker threads)
static void processChunk(const std::string& name,
                         const std::string& seq,
                         int chunkStart, int chunkStop,
                         std::vector<BamFile*>& bamFiles,
                         std::vector<GenomeRegion>& results,
                         std::mutex& resultsMutex,
                         std::atomic<int>& completedChunks,
                         int totalChunks) {
    // Each thread opens its own BAM file handles (htslib is NOT thread-safe)
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

    // Process each BAM file
    for (auto* bam : threadBams) {
        bam->process(region);
    }

    // Store result in thread-safe manner
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

    // Clean up thread-local BAM handles
    for (auto* bam : threadBams) {
        bam->close();
        delete bam;
    }
}

GenomeRegion::GenomeRegion(const std::string& name, int start, int stop,
                           const std::string& bases, double minDepth)
    : name(name), start(start), stop(stop), contigBases(bases),
      originalBases(bases), minDepth(minDepth) {
    pileUps.resize(stop - start);
}

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
    return baseAt(pos);
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

const PileUp & GenomeRegion::pileUpRegion(int index) const {
    return pileUps[index];
}

int GenomeRegion::locus(int index) const {
    return start + index;
}

// GenomeFile implementation
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
            // Parse header - take first word after >
            std::istringstream iss(line.substr(1));
            iss >> currentName;
        } else {
            // Convert to uppercase
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
            // Find matching contig
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
        // Single-threaded mode (original behavior)
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
                    if (bam) {
                        bam->process(region);
                    }
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

    // Collect all chunks across all contigs
    struct ChunkTask {
        std::string name;
        std::string seq;  // We store the full seq per contig, substr in worker
        int chunkStart;
        int chunkStop;
        int contigLength;
    };

    std::vector<ChunkTask> tasks;
    for (const auto& contig : contigs) {
        const std::string& name = contig.first;
        const std::string& seq = contig.second;
        int length = static_cast<int>(seq.size());

        std::cout << "Processing " << name << " (" << length << " bp)" << std::endl;

        for (int chunkStart = 0; chunkStart < length; chunkStart += Pilon::chunkSize) {
            int chunkStop = std::min(chunkStart + Pilon::chunkSize, length);
            tasks.push_back({name, seq, chunkStart, chunkStop, length});
        }
    }

    int totalChunks = static_cast<int>(tasks.size());
    std::cout << "Total chunks: " << totalChunks << std::endl;

    // Thread-safe results storage
    std::mutex resultsMutex;
    std::atomic<int> completedChunks(0);

    // Pre-allocate results vector
    processedRegions_.reserve(totalChunks);

    // Launch worker threads
    std::vector<std::thread> workers;
    workers.reserve(numThreads);

    // Simple work-stealing queue with mutex
    std::mutex taskMutex;
    std::queue<size_t> taskQueue;
    for (size_t i = 0; i < tasks.size(); i++) {
        taskQueue.push(i);
    }

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

    for (int i = 0; i < numThreads; i++) {
        workers.emplace_back(workerFunc);
    }

    for (auto& t : workers) {
        t.join();
    }

    // Sort results by genomic position (chunks may have completed out of order)
    std::sort(processedRegions_.begin(), processedRegions_.end(),
              [](const GenomeRegion& a, const GenomeRegion& b) {
                  if (a.name != b.name) return a.name < b.name;
                  return a.start < b.start;
              });

    std::cout << "Multi-threaded processing complete: " << totalChunks << " chunks" << std::endl;
}

} // namespace pilon
