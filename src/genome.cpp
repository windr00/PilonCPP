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

namespace pilon {

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

    for (const auto& contig : contigs) {
        const std::string& name = contig.first;
        const std::string& seq = contig.second;
        int length = static_cast<int>(seq.size());

        std::cout << "Processing " << name << " (" << length << " bp)" << std::endl;

        // Process in chunks
        for (int chunkStart = 0; chunkStart < length; chunkStart += Pilon::chunkSize) {
            int chunkStop = std::min(chunkStart + Pilon::chunkSize, length);

            GenomeRegion region(name, chunkStart, chunkStop,
                               seq.substr(chunkStart, chunkStop - chunkStart),
                               Pilon::minDepth);

            // Process each BAM file
            for (auto* bam : bamFiles) {
                if (bam) {
                    bam->process(region);
                }
            }

            // Store processed region for later VCF output
            processedRegions_.push_back(std::move(region));

            if (Pilon::verbose) {
                std::cout << "  Chunk " << chunkStart << "-" << chunkStop << " done" << std::endl;
            }
        }
    }
}

} // namespace pilon
