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
const std::unordered_set<std::string> Pilon::fixChoices = {"all", "snps", "indels", "gaps", "local", "bases", "none"};
const std::unordered_set<std::string> Pilon::experimentalFixChoices = {"amb", "breaks", "circles", "novel", "scaffolds"};

bool Pilon::fixSnps = true;
bool Pilon::fixIndels = true;
bool Pilon::fixGaps = true;
bool Pilon::fixLocal = true;
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

int Pilon::chunkSize = 10000000;
uint8_t Pilon::defaultQual = 10;
bool Pilon::diploid = false;
bool Pilon::duplicates = false;
bool Pilon::dumpReads = false;
std::unordered_set<std::string> Pilon::fixList;
int Pilon::flank = 10;
int Pilon::gapMargin = 100000;
bool Pilon::iupac = false;
int Pilon::minMinDepth = 5;
int Pilon::minGap = 10;
double Pilon::minDepth = 0.1;
int Pilon::minQual = 0;
int Pilon::minMq = 0;
bool Pilon::multiClosure = false;
bool Pilon::nonPf = false;
bool Pilon::oldIndel = false;
bool Pilon::longread = false;
bool Pilon::pacbio = false;
bool Pilon::nanopore = false;
bool Pilon::strays = true;
bool Pilon::trSafe = true;

std::vector<std::string> Pilon::novelContigs;

// Threading
int Pilon::threads = 1;

// For logging
std::vector<std::string> Pilon::commandArgs;

void Pilon::parseFixList(const std::string& fix) {
    auto fixes = Utils::split(fix, ',');
    for (const auto& f : fixes) {
        std::string trimmed = Utils::trim(f);
        if (fixChoices.count(trimmed)) {
            fixList.insert(trimmed);
            if (trimmed == "all" || trimmed == "everything") {
                fixSnps = true;
                fixIndels = true;
                fixGaps = true;
                fixLocal = true;
            } else if (trimmed == "bases") {
                fixSnps = true;
                fixIndels = true;
            } else if (trimmed == "none") {
                fixSnps = false;
                fixIndels = false;
                fixGaps = false;
                fixLocal = false;
            } else if (trimmed == "snps") fixSnps = true;
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
    std::cerr << std::endl;
    std::cerr << "PilonCpp version 1.0.0" << std::endl;
    std::cerr << std::endl;
    std::cerr << "    Usage: piloncpp --genome genome.fasta [--frags frags.bam] [--jumps jumps.bam] [--unpaired unpaired.bam]" << std::endl;
    std::cerr << "                 [...other options...]" << std::endl;
    std::cerr << "           piloncpp --help for option details" << std::endl;
}

void Pilon::printHelp() {
    printUsage();
    std::cerr << R"(
         INPUTS:
           --genome genome.fasta (or --input)
              The input genome we are trying to improve, which must be the reference used
              for the bam alignments.  At least one of --frags or --jumps must also be given.
           --frags frags.bam
              A bam file consisting of fragment paired-end alignments, aligned to the --genome
              argument using bwa or bowtie2.  This argument may be specifed more than once.
           --jumps jumps.bam
              A bam file consisting of jump (mate pair) paired-end alignments, aligned to the
              --genome argument using bwa or bowtie2.  This argument may be specifed more than once.
           --unpaired unpaired.bam
              A bam file consisting of unpaired alignments, aligned to the --genome argument
              using bwa or bowtie2.  This argument may be specifed more than once.
           --bam any.bam
              A bam file of unknown type; Pilon will scan it and attempt to classify it as one
              of the above bam types.
           --nanopore ont.bam
              A bam file containing Oxford Nanopore read alignments. Experimental.
           --pacbio pb.bam
              A bam file containing Pacific Biosciences read alignments. Experimental.
         OUTPUTS:
           --output prefix (or --prefix)
              Prefix for output files
           --outdir directory
              Use this directory for all output files.
           --changes
              If specified, a file listing changes in the <output>.fasta will be generated.
           --vcf
              If specified, a vcf file will be generated
           --vcfqe (or --vcfQE)
               If specified, the VCF will contain a QE (quality-weighted evidence) field rather
               than the default QP (quality-weighted percentage of evidence) field.
           --tracks
               This options will cause many track files (*.bed, *.wig) suitable for viewing in
               a genome browser to be written.
         CONTROL:
           --variant
              Sets up heuristics for variant calling, as opposed to assembly improvement;
              equivalent to "--vcf --fix all,breaks".
           --chunksize (or --chunk-size)
              Input FASTA elements larger than this will be processed in smaller pieces not to
              exceed this size (default 10000000).
           --diploid
              Sample is from diploid organism; will eventually affect calling of heterozygous SNPs
           --fix fixlist
              A comma-separated list of categories of issues to try to fix:
                "snps": try to fix individual base errors;
                "indels": try to fix small indels;
                "gaps": try to fill gaps;
                "local": try to detect and fix local misassemblies;
                "all": all of the above (default);
                "bases": shorthand for "snps" and "indels" (for back compatibility);
                "none": none of the above; new fasta file will not be written.
              The following are experimental fix types:
                "amb": fix ambiguous bases in fasta output (to most likely alternative);
                "breaks": allow local reassembly to open new gaps (with "local");
                "circles": try to close circlar elements when used with long corrected reads;
                "novel": assemble novel sequence from unaligned non-jump reads.
           --dumpreads
              Dump reads for local re-assemblies.
           --duplicates
              Use reads marked as duplicates in the input BAMs (ignored by default).
           --iupac
              Output IUPAC ambiguous base codes in the output FASTA file when appropriate.
           --nonpf
              Use reads which failed sequencer quality filtering (ignored by default).
           --nostrays
              Skip making a pass through the input BAM files to identify stray pairs.
           --targets targetlist
              Only process the specified target(s).  Targets are comma-separated, and each target
              is a fasta element name optionally followed by a base range.
              Example: "scaffold00001,scaffold00002:10000-20000" would result in processing all of
              scaffold00001 and coordinates 10000-20000 of scaffold00002.
              If "targetlist" is the name of a file, each line will be treated as a target
              specification.
           --threads
              Number of parallel threads to use (default: 1).
              Note: this option is C++ specific; the Scala version does not support multithreading.
           --verbose
              More verbose output.
           --debug
              Debugging output (implies verbose).
           --version
              Print version string and exit.
           --help, -h
              Show this help message.
         HEURISTICS:
           --defaultqual qual
              Assumes bases are of this quality if quals are no present in input BAMs (default 10).
           --flank nbases
              Controls how much of the well-aligned reads will be used; this many bases at each
              end of the good reads will be ignored (default 10).
           --gapmargin (or --gap-margin)
              Closed gaps must be within this number of bases of true size to be closed (100000)
           --kmer (or --K)
              Kmer size used by internal assembler (default 47).
           --mindepth (or --min-depth)
              Variants (snps and indels) will only be called if there is coverage of good pairs
              at this depth or more; if this value is >= 1, it is an absolute depth, if it is a
              fraction < 1, then minimum depth is computed by multiplying this value by the mean
              coverage for the region, with a minumum value of 5 (default 0.1: min depth to call
              is 10% of mean coverage or 5, whichever is greater).
           --mingap (or --min-gap)
              Minimum size for unclosed gaps (default 10)
           --minmq (or --min-mq)
              Minimum alignment mapping quality for a read to count in pileups (default 0)
           --minqual (or --min-qual)
              Minimum base quality to consider for pileups (default 0)
)";
}

void Pilon::parseOptions(int argc, char* argv[]) {
    // Store command line arguments
    for (int i = 0; i < argc; i++) {
        commandArgs.push_back(argv[i]);
    }

    static struct option longOptions[] = {
        {"input",         required_argument, 0, 'i'},
        {"genome",        required_argument, 0, 'i'},
        {"output",        required_argument, 0, 'o'},
        {"prefix",        required_argument, 0, 'o'},
        {"frags",         required_argument, 0, 'f'},
        {"jumps",         required_argument, 0, 'j'},
        {"unpaired",      required_argument, 0, 'U'},
        {"bam",           required_argument, 0, 'b'},
        {"nanopore",      required_argument, 0, 'N'},
        {"pacbio",        required_argument, 0, 'P'},
        {"fix",           required_argument, 0, 'x'},
        {"changes",       no_argument,       0, 'c'},
        {"tracks",        no_argument,       0, 't'},
        {"vcf",           no_argument,       0, 'v'},
        {"vcfqe",         no_argument,       0, 'q'},
        {"vcfQE",         no_argument,       0, 'q'},
        {"diploid",       no_argument,       0, 'd'},
        {"haploid",       no_argument,       0, 'p'},
        {"duplicates",    no_argument,       0, 'u'},
        {"no-duplicates", no_argument,       0, 'n'},
        {"nostrays",      no_argument,       0, 'S'},
        {"dumpreads",     no_argument,       0, 'R'},
        {"iupac",         no_argument,       0, 'I'},
        {"nonpf",         no_argument,       0, 'K'},
        {"variant",       no_argument,       0, 'W'},
        {"min-depth",     required_argument, 0, 'm'},
        {"mindepth",      required_argument, 0, 'm'},
        {"min-qual",      required_argument, 0, 'l'},
        {"minqual",       required_argument, 0, 'l'},
        {"min-mq",        required_argument, 0, 'M'},
        {"minmq",         required_argument, 0, 'M'},
        {"min-min-depth", required_argument, 0, 'D'},
        {"min-min-depth", required_argument, 0, 'D'},
        {"chunk-size",    required_argument, 0, 's'},
        {"chunksize",     required_argument, 0, 's'},
        {"flank",         required_argument, 0, 'F'},
        {"gap-margin",    required_argument, 0, 'G'},
        {"gapmargin",     required_argument, 0, 'G'},
        {"min-gap",       required_argument, 0, 'g'},
        {"mingap",        required_argument, 0, 'g'},
        {"defaultqual",   required_argument, 0, 'Q'},
        {"kmer",          required_argument, 0, 2},
        {"K",             required_argument, 0, 2},
        {"outdir",        required_argument, 0, 'O'},
        {"targets",       required_argument, 0, 'T'},
        {"verbose",       no_argument,       0, 'V'},
        {"debug",         no_argument,       0, 'Z'},
        {"version",       no_argument,       0, 3},
        {"threads",       required_argument, 0, 1},
        {"help",          no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    int optionIndex = 0;

    while ((opt = getopt_long(argc, argv, "i:o:f:j:b:U:N:P:x:ctvqdpunSRIm:l:M:D:s:F:G:g:Q:O:T:VZh", longOptions, &optionIndex)) != -1) {
        switch (opt) {
            case 'i': genomePath = optarg; break;
            case 'o': prefix = optarg; break;
            case 'f':
            case 'j':
            case 'U':
            case 'N':
            case 'P':
            case 'b': {
                std::string type;
                if (opt == 'f') type = "frags";
                else if (opt == 'j') type = "jumps";
                else if (opt == 'U') type = "unpaired";
                else if (opt == 'N') { type = "nanopore"; nanopore = true; longread = true; }
                else if (opt == 'P') { type = "pacbio"; pacbio = true; longread = true; }
                else type = "bam";
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
            case 'S': strays = false; break;
            case 'R': dumpReads = true; break;
            case 'I': iupac = true; break;
            case 'K': nonPf = true; break;
            case 'W': {
                // --variant: equvialent to --vcf --fix all,breaks
                vcf = true;
                parseFixList("all,breaks");
                break;
            }
            case 'm': minDepth = std::stod(optarg); break;
            case 'l': minQual = std::stoi(optarg); break;
            case 'M': minMq = std::stoi(optarg); break;
            case 'D': minMinDepth = std::stoi(optarg); break;
            case 's': chunkSize = std::stoi(optarg); break;
            case 'F': flank = std::stoi(optarg); break;
            case 'G': gapMargin = std::stoi(optarg); break;
            case 'g': minGap = std::stoi(optarg); break;
            case 'Q': defaultQual = static_cast<uint8_t>(std::stoi(optarg)); break;
            case 2:   /* --kmer / --K: stored but not used yet */ break;
            case 'O': outdir = optarg; break;
            case 'T': targets = optarg; break;
            case 'V': verbose = true; break;
            case 'Z': debug = true; verbose = true; break;
            case 1:   threads = std::stoi(optarg); break;
            case 3:   // --version
                std::cout << "PilonCpp version 1.0.0" << std::endl;
                exit(0);
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
        // Scala default: all fixes enabled
        fixList = fixChoices;
    }

    // Scala: strays is only active when fixing gaps/local/scaffolds
    strays = strays && (fixGaps || fixLocal || fixScaffolds);
}

} // namespace pilon
