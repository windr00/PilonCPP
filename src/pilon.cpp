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
#include "utils.h"
#include "bamfile.h"
#include <iostream>
#include <cstring>
#include <getopt.h>

namespace pilon {

// Static member initialization
const std::unordered_set<std::string> Pilon::fixChoices = {"snps", "indels", "gaps", "local"};
const std::unordered_set<std::string> Pilon::experimentalFixChoices = {"amb", "breaks", "circles", "novel", "scaffolds"};

bool Pilon::fixSnps = false;
bool Pilon::fixIndels = false;
bool Pilon::fixGaps = false;
bool Pilon::fixLocal = false;
bool Pilon::fixAmb = false;
bool Pilon::fixBreaks = false;
bool Pilon::fixCircles = false;
bool Pilon::fixNovel = false;
bool Pilon::fixScaffolds = false;

std::vector<BamFile*> Pilon::bamFiles;
std::string Pilon::targets;
std::string Pilon::genomePath;
std::string Pilon::prefix;
std::string Pilon::outdir;
bool Pilon::changes = false;
bool Pilon::tracks = false;
bool Pilon::verbose = false;
bool Pilon::vcf = false;
bool Pilon::vcfQE = false;
bool Pilon::debug = false;

int Pilon::chunkSize = 100000;
uint8_t Pilon::defaultQual = 30;
bool Pilon::diploid = true;
bool Pilon::duplicates = true;
bool Pilon::dumpReads = false;
std::unordered_set<std::string> Pilon::fixList;
int Pilon::flank = 1000;
int Pilon::gapMargin = 100;
bool Pilon::iupac = false;
int Pilon::minMinDepth = 10;
int Pilon::minGap = 10;
double Pilon::minDepth = 10.0;
int Pilon::minQual = 13;
int Pilon::minMq = 5;
bool Pilon::multiClosure = false;
bool Pilon::nonPf = false;
bool Pilon::oldIndel = false;
bool Pilon::longread = false;
bool Pilon::pacbio = false;
bool Pilon::nanopore = false;
bool Pilon::strays = false;
bool Pilon::trSafe = false;

std::vector<std::string> Pilon::novelContigs;
std::vector<std::string> Pilon::commandArgs;

void Pilon::parseFixList(const std::string& fix) {
    auto fixes = Utils::split(fix, ',');
    for (const auto& f : fixes) {
        std::string trimmed = Utils::trim(f);
        if (fixChoices.count(trimmed)) {
            fixList.insert(trimmed);
            if (trimmed == "snps") fixSnps = true;
            else if (trimmed == "indels") fixIndels = true;
            else if (trimmed == "gaps") fixGaps = true;
            else if (trimmed == "local") fixLocal = true;
        } else if (experimentalFixChoices.count(trimmed)) {
            fixList.insert(trimmed);
            if (trimmed == "amb") fixAmb = true;
            else if (trimmed == "breaks") fixBreaks = true;
            else if (trimmed == "circles") fixCircles = true;
            else if (trimmed == "novel") fixNovel = true;
            else if (trimmed == "scaffolds") fixScaffolds = true;
        } else {
            std::cerr << "Unknown fix: " << trimmed << std::endl;
        }
    }
}

std::string Pilon::outputFile(const std::string& name) {
    if (!outdir.empty()) return outdir + "/" + name;
    if (!prefix.empty()) return prefix + "." + name;
    return name;
}

void Pilon::printUsage() {
    std::cerr << "Usage: piloncpp [options] --input <genome> --output <output>" << std::endl;
    std::cerr << "       --frags <bam1> [--frags <bam2> ...]" << std::endl;
}

void Pilon::printHelp() {
    printUsage();
    std::cerr << R"(
Options:
  --input, -i       Input genome (FASTA)
  --output, -o      Output genome (FASTA)
  --frags           Fragment BAM files
  --jumps           Jumping library BAM files
  --bam             BAM files (auto-detect type)
  --fix             Comma-separated list of fixes: snps,indels,gaps,local
                    Experimental: amb,breaks,circles,novel,scaffolds
  --changes         Output only changes
  --tracks          Output tracks for IGV
  --vcf             Output VCF file
  --vcfQE           Output VCF with quality estimates
  --diploid         Assume diploid genome (default)
  --haploid         Assume haploid genome
  --duplicates      Include duplicate reads (default)
  --no-duplicates   Exclude duplicate reads
  --min-depth       Minimum depth for calling (default: 10)
  --min-qual        Minimum base quality (default: 13)
  --min-mq          Minimum mapping quality (default: 5)
  --min-min-depth   Minimum minimum depth (default: 10)
  --chunk-size      Chunk size for processing (default: 100000)
  --flank           Flank size for gap filling (default: 1000)
  --gap-margin      Gap margin (default: 100)
  --min-gap         Minimum gap size to fill (default: 10)
  --outdir          Output directory
  --prefix          Output prefix
  --targets         Target regions file
  --verbose         Verbose output
  --debug           Debug output
  --help, -h        Show this help message
)";
}

void Pilon::parseOptions(int argc, char* argv[]) {
    // Store command line arguments
    for (int i = 0; i < argc; i++) {
        commandArgs.push_back(argv[i]);
    }

    static struct option longOptions[] = {
        {"input",         required_argument, 0, 'i'},
        {"output",        required_argument, 0, 'o'},
        {"frags",         required_argument, 0, 'f'},
        {"jumps",         required_argument, 0, 'j'},
        {"bam",           required_argument, 0, 'b'},
        {"fix",           required_argument, 0, 'x'},
        {"changes",       no_argument,       0, 'c'},
        {"tracks",        no_argument,       0, 't'},
        {"vcf",           no_argument,       0, 'v'},
        {"vcfQE",         no_argument,       0, 'q'},
        {"diploid",       no_argument,       0, 'd'},
        {"haploid",       no_argument,       0, 'p'},
        {"duplicates",    no_argument,       0, 'u'},
        {"no-duplicates", no_argument,       0, 'n'},
        {"min-depth",     required_argument, 0, 'm'},
        {"min-qual",      required_argument, 0, 'l'},
        {"min-mq",        required_argument, 0, 'M'},
        {"min-min-depth", required_argument, 0, 'D'},
        {"chunk-size",    required_argument, 0, 's'},
        {"flank",         required_argument, 0, 'F'},
        {"gap-margin",    required_argument, 0, 'G'},
        {"min-gap",       required_argument, 0, 'g'},
        {"outdir",        required_argument, 0, 'O'},
        {"prefix",        required_argument, 0, 'P'},
        {"targets",       required_argument, 0, 'T'},
        {"verbose",       no_argument,       0, 'V'},
        {"debug",         no_argument,       0, 'Z'},
        {"help",          no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    int optionIndex = 0;

    while ((opt = getopt_long(argc, argv, "i:o:f:j:b:x:ctvqdpnul:M:D:s:F:G:g:O:P:T:VZh", longOptions, &optionIndex)) != -1) {
        switch (opt) {
            case 'i': genomePath = optarg; break;
            case 'o': prefix = optarg; break;
            case 'f':
            case 'j':
            case 'b': {
                std::string type = (opt == 'f') ? "frags" : (opt == 'j') ? "jumps" : "bam";
                bamFiles.push_back(new BamFile(optarg, type));
                break;
            }
            case 'x': parseFixList(optarg); break;
            case 'c': changes = true; break;
            case 't': tracks = true; break;
            case 'v': vcf = true; break;
            case 'q': vcfQE = true; break;
            case 'd': diploid = true; break;
            case 'p': diploid = false; break;
            case 'u': duplicates = true; break;
            case 'n': duplicates = false; break;
            case 'm': minDepth = std::stod(optarg); break;
            case 'l': minQual = std::stoi(optarg); break;
            case 'M': minMq = std::stoi(optarg); break;
            case 'D': minMinDepth = std::stoi(optarg); break;
            case 's': chunkSize = std::stoi(optarg); break;
            case 'F': flank = std::stoi(optarg); break;
            case 'G': gapMargin = std::stoi(optarg); break;
            case 'g': minGap = std::stoi(optarg); break;
            case 'O': outdir = optarg; break;
            case 'P': prefix = optarg; break;
            case 'T': targets = optarg; break;
            case 'V': verbose = true; break;
            case 'Z': debug = true; break;
            case 'h':
                printHelp();
                exit(0);
            default:
                printUsage();
                exit(1);
        }
    }

    // Validate required options
    if (genomePath.empty()) {
        std::cerr << "Error: --input is required" << std::endl;
        printUsage();
        exit(1);
    }
    if (prefix.empty()) {
        std::cerr << "Error: --output is required" << std::endl;
        printUsage();
        exit(1);
    }
    if (bamFiles.empty()) {
        std::cerr << "Error: at least one BAM file is required" << std::endl;
        printUsage();
        exit(1);
    }
    if (fixList.empty()) {
        std::cerr << "Error: --fix is required (e.g., --fix snps,indels)" << std::endl;
        printUsage();
        exit(1);
    }
}

} // namespace pilon
