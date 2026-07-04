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

#ifndef PILON_PILON_H
#define PILON_PILON_H

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_set>

namespace pilon {

class BamFile;

class Pilon {
public:
    // Fix choices
    static const std::unordered_set<std::string> fixChoices;
    static const std::unordered_set<std::string> experimentalFixChoices;

    // Fix flags
    static bool fixSnps;
    static bool fixIndels;
    static bool fixGaps;
    static bool fixLocal;
    static bool fixAmb;
    static bool fixBreaks;
    static bool fixCircles;
    static bool fixNovel;
    static bool fixScaffolds;

    // Input parameters
    static std::vector<BamFile*> bamFiles;
    static std::string targets;
    static std::string genomePath;

    // Output parameters
    static std::string prefix;
    static std::string outdir;
    static bool changes;
    static bool tracks;
    static bool verbose;
    static bool vcf;
    static bool vcfQE;
    static bool debug;

    // Heuristics and control parameters
    static int chunkSize;
    static uint8_t defaultQual;
    static bool diploid;
    static bool duplicates;
    static bool dumpReads;
    static std::unordered_set<std::string> fixList;
    static int flank;
    static int gapMargin;
    static bool iupac;
    static int minMinDepth;
    static int minGap;
    static double minDepth;
    static int minQual;
    static int minMq;
    static bool multiClosure;
    static bool nonPf;
    static bool oldIndel;
    static bool longread;
    static bool pacbio;
    static bool nanopore;
    static bool strays;
    static bool trSafe;

    // Global computed data
    static std::vector<std::string> novelContigs;

    // Threading
    static int threads;

    // For logging
    static std::vector<std::string> commandArgs;

    static void parseOptions(int argc, char* argv[]);
    static void parseFixList(const std::string& fix);
    static std::string outputFile(const std::string& name);
    static void printUsage();
    static void printHelp();
};

} // namespace pilon

#endif // PILON_PILON_H
