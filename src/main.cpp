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
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace pilon;

int main(int argc, char* argv[]) {
    // Parse command line options
    Pilon::parseOptions(argc, argv);

    std::cout << "PilonCpp - Genome assembly polishing tool" << std::endl;
    std::cout << "=========================================" << std::endl;
    std::cout << "Input genome: " << Pilon::genomePath << std::endl;
    std::cout << "Output prefix: " << Pilon::prefix << std::endl;
    std::cout << "BAM files: " << Pilon::bamFiles.size() << std::endl;
    std::cout << "Fixes: ";
    for (const auto& fix : Pilon::fixList) {
        std::cout << fix << " ";
    }
    std::cout << std::endl;

    // Load genome
    std::cout << "\nLoading genome..." << std::endl;
    GenomeFile genome(Pilon::genomePath, Pilon::targets);
    
    auto contigs = genome.getContigs();
    std::cout << "Loaded " << contigs.size() << " contigs" << std::endl;

    // Open BAM files
    std::cout << "\nOpening BAM files..." << std::endl;
    for (auto* bam : Pilon::bamFiles) {
        if (bam->open()) {
            std::cout << "  Opened: " << bam->path() << " (" << bam->bamType() << ")" << std::endl;
        } else {
            std::cerr << "  Failed to open: " << bam->path() << std::endl;
        }
    }

    // Scan BAM files for statistics
    std::cout << "\nScanning BAM files..." << std::endl;
    auto seqNames = genome.getContigs();
    std::unordered_set<std::string> seqsOfInterest;
    for (const auto& c : seqNames) {
        seqsOfInterest.insert(c.first);
    }
    
    for (auto* bam : Pilon::bamFiles) {
        if (bam) {
            bam->scan(seqsOfInterest);
        }
    }

    // Process genome regions
    std::cout << "\nProcessing genome regions..." << std::endl;
    genome.processRegions(Pilon::bamFiles);

    // Output VCF if requested
    if (Pilon::vcf) {
        std::string vcfPath = Pilon::outputFile("variants.vcf");
        std::cout << "\nWriting VCF: " << vcfPath << std::endl;
        
        auto contigSizes = genome.getContigSizes();
        Vcf vcf(vcfPath, contigSizes);
        
        // Write variant records from processed regions
        const auto& regions = genome.getProcessedRegions();
        int totalVariants = 0;
        
        for (const auto& region : regions) {
            for (int i = 0; i < region.size(); i++) {
                const PileUp& pu = region.pileUpRegion(i);
                auto bc = pu.baseCall();
                
                // Skip if no call or depth too low
                if (!bc.called() || pu.depth() < Pilon::minDepth) continue;
                
                char refBase = region.refBase(i);
                
                // Write SNP variants
                if (Pilon::fixSnps && bc.called() && !bc.isInsertion() && !bc.isDeletion()) {
                    if (!bc.baseMatch(refBase)) {
                        vcf.writeRecord(region, i, false, false);
                        totalVariants++;
                    }
                }
                
                // Write indel variants
                if (Pilon::fixIndels && (bc.isInsertion() || bc.isDeletion())) {
                    vcf.writeRecord(region, i, false, true);
                    totalVariants++;
                }
            }
        }
        
        vcf.close();
        std::cout << "  Written " << totalVariants << " variant records" << std::endl;
    }

    // Output polished genome
    std::string outputPath = Pilon::outputFile("fasta");
    std::cout << "\nWriting polished genome: " << outputPath << std::endl;
    
    std::ofstream outFile(outputPath);
    if (!outFile.is_open()) {
        std::cerr << "Error: Cannot open output file: " << outputPath << std::endl;
        return 1;
    }

    for (const auto& contig : contigs) {
        outFile << ">" << contig.first << std::endl;
        // Write sequence in 60-character lines
        const std::string& seq = contig.second;
        for (size_t i = 0; i < seq.length(); i += 60) {
            outFile << seq.substr(i, 60) << std::endl;
        }
    }
    outFile.close();

    // Clean up
    for (auto* bam : Pilon::bamFiles) {
        delete bam;
    }
    Pilon::bamFiles.clear();

    std::cout << "\nDone!" << std::endl;
    return 0;
}
