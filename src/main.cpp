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

#include "pilon.h"
#include "genome.h"
#include "vcf.h"
#include "gapfiller.h"
#include "bamfile.h"
#include "utils.h"
#include "bases.h"
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace pilon;

int main(int argc, char* argv[]) {
    Pilon::parseOptions(argc, argv);

    if (Pilon::outdir != "") {
        Utils::mkdirs(Pilon::outdir);
    }

    // Scala: strays &= fixGaps || fixLocal || fixScaffolds (done in parseOptions)

    std::cout << "PilonCpp - Genome assembly polishing tool" << std::endl;
    std::cout << "Genome: " << Pilon::genomePath << std::endl;
    std::cout << "Fixing " << Pilon::fixList.size() << " categories" << std::endl;

    GenomeFile genome(Pilon::genomePath, Pilon::targets);

    // Open BAM files
    for (auto* bam : Pilon::bamFiles) {
        if (!bam->open()) {
            std::cerr << "Failed to open BAM: " << bam->path() << std::endl;
        }
    }

    // Scan BAM files (matching Scala: scan before processing)
    if (Pilon::strays || Pilon::fixCircles) {
        std::cout << "Scanning BAMs" << std::endl;
        auto contigs = genome.getContigs();
        std::unordered_set<std::string> seqsOfInterest;
        for (const auto& c : contigs) seqsOfInterest.insert(c.first);
        for (auto* bam : Pilon::bamFiles) {
            if (bam) bam->scan(seqsOfInterest);
        }
    }

    // fixNovel: assemble novel contigs before region processing
    if (Pilon::fixNovel) {
        std::cout << "Assembling novel sequence" << std::endl;
        Assembler genomeGraph(1);
        for (const auto& contig : genome.getContigs()) {
            genomeGraph.addGraphSeq(contig.second);
        }
        Assembler assembler;
        for (auto* bam : Pilon::bamFiles) {
            if (bam && bam->bamType() != "jumps") {
                auto reads = bam->getUnalignedReads();
                assembler.addReads(reads);
            }
        }
        auto novelContigs = assembler.novel(genomeGraph);
        Pilon::novelContigs = novelContigs;

        int totalBases = 0;
        for (const auto& nc : novelContigs) totalBases += nc.size();
        std::cout << "Assembled " << novelContigs.size() << " novel contigs containing "
                  << totalBases << " bases" << std::endl;
    }

    // Process genome regions
    std::cout << "Input genome size: " << genome.getContigSizes().size() << " contigs" << std::endl;
    genome.processRegions(Pilon::bamFiles);

    // Set up VCF if requested (matching Scala)
    Vcf* vcf = nullptr;
    if (Pilon::vcf) {
        auto contigSizes = genome.getContigSizes();
        std::string vcfPath = Pilon::outputFile(".vcf");
        std::cout << "Writing VCF to " << vcfPath << std::endl;
        vcf = new Vcf(vcfPath, contigSizes);
    }

    // Set up changes file if requested
    std::ofstream* changesFile = nullptr;
    if (Pilon::changes) {
        std::string changesPath = Pilon::outputFile(".changes");
        changesFile = new std::ofstream(changesPath);
        std::cout << "Writing changes to " << changesPath << std::endl;
    }

    // Write polished FASTA (matching Scala)
    bool writeFasta = Pilon::fixSnps || Pilon::fixIndels || Pilon::fixGaps ||
                      Pilon::fixLocal || Pilon::fixNovel;
    std::ofstream* fastaFile = nullptr;
    if (writeFasta) {
        std::string fastaPath = Pilon::outputFile(".fasta");
        std::cout << "Writing polished genome to " << fastaPath << std::endl;
        fastaFile = new std::ofstream(fastaPath);
    }

    auto& regions = genome.getProcessedRegions();
    std::string currentName;
    std::string currentSeq;

    for (const auto& region : regions) {
        if (region.name != currentName) {
            // Write previous contig
            if (!currentName.empty() && fastaFile && fastaFile->is_open()) {
                std::string newName = currentName;
                if (currentName.find('|') == std::string::npos) {
                    newName += "_pilon";
                } else if (currentName.back() != '|') {
                    newName += "|pilon";
                } else {
                    newName += "pilon";
                }
                (*fastaFile) << ">" << newName << "\n";
                for (size_t i = 0; i < currentSeq.length(); i += 80) {
                    (*fastaFile) << currentSeq.substr(i, 80) << "\n";
                }
            }
            currentName = region.name;
            currentSeq.clear();
        }
        currentSeq += region.bases;

        // Write VCF records for this region
        if (vcf) {
            // Write fix records for local reassembly (matching Scala writeFixRecord)
            // Scala only writes for big fixes (reassembly), not pileup-based small indels
            for (const auto& fix : region.fixes) {
                const std::string& ref = std::get<1>(fix);
                const std::string& patch = std::get<2>(fix);
                // Big fix: involves large change or contains N (gap fill/break repair)
                if (ref.length() >= 10 || patch.length() >= 10 ||
                    ref.find('N') != std::string::npos || patch.find('N') != std::string::npos) {
                    vcf->writeFixRecord(region, fix);
                }
            }
            // Write per-position records
            for (int i = 0; i < region.size(); i++) {
                vcf->writeRecord(region, i, false, true);
            }
        }

        // Write changes
        if (changesFile && changesFile->is_open()) {
            int delta = 0;
            for (const auto& fix : region.fixes) {
                int loc = std::get<0>(fix) + delta;
                const std::string& from = std::get<1>(fix);
                const std::string& to = std::get<2>(fix);
                std::string fromDisplay = from.empty() ? "." : from;
                std::string toDisplay = to.empty() ? "." : to;
                (*changesFile) << region.name << "\t" << loc << "\t" << fromDisplay << "\t" << toDisplay << "\n";
                delta += static_cast<int>(to.length()) - static_cast<int>(from.length());
            }
        }
    }

    // Write last contig
    if (!currentName.empty() && fastaFile && fastaFile->is_open()) {
        std::string newName = currentName;
        if (currentName.find('|') == std::string::npos) {
            newName += "_pilon";
        } else if (currentName.back() != '|') {
            newName += "|pilon";
        }
        (*fastaFile) << ">" << newName << "\n";
        for (size_t i = 0; i < currentSeq.length(); i += 80) {
            (*fastaFile) << currentSeq.substr(i, 80) << "\n";
        }
    }

    // Write novel contigs to FASTA (matching Scala)
    if (Pilon::fixNovel && fastaFile && fastaFile->is_open()) {
        const auto& novelContigs = Pilon::novelContigs;
        for (size_t n = 0; n < novelContigs.size(); n++) {
            char header[64];
            snprintf(header, sizeof(header), "pilon_novel_%03zu", n + 1);
            std::cout << "Appending " << header << " length " << novelContigs[n].size() << std::endl;
            (*fastaFile) << ">" << header << "\n";
            for (size_t i = 0; i < novelContigs[n].size(); i += 80) {
                (*fastaFile) << novelContigs[n].substr(i, 80) << "\n";
            }
        }
    }

    if (fastaFile) { fastaFile->close(); delete fastaFile; }
    if (vcf) { vcf->close(); delete vcf; }
    if (changesFile) { changesFile->close(); delete changesFile; }

    // Output IGV tracks
    if (Pilon::tracks) {
        std::string trackPrefix = Pilon::outputFile("");
        std::cout << "Writing tracks" << std::endl;
        for (auto& region : regions) {
            region.writeTracks(trackPrefix);
        }
    }

    // Cleanup
    for (auto* bam : Pilon::bamFiles) {
        delete bam;
    }
    Pilon::bamFiles.clear();

    std::cout << "Done!" << std::endl;
    return 0;
}
