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

#include "vcf.h"
#include "genome.h"
#include "pilon.h"
#include <iostream>

namespace pilon {

Vcf::Vcf(const std::string& filePath,
         const std::vector<std::pair<std::string, int>>& contigsWithSizes)
    : filePath_(filePath) {
    writer_ = fopen(filePath.c_str(), "w");
    if (!writer_) {
        std::cerr << "Error: Cannot open VCF file: " << filePath << std::endl;
        exit(1);
    }
    writeHeader(contigsWithSizes);
}

Vcf::~Vcf() {
    close();
}

void Vcf::writeHeader(const std::vector<std::pair<std::string, int>>& contigsWithSizes) {
    fprintf(writer_, "##fileformat=VCFv4.2\n");
    fprintf(writer_, "##source=PilonCpp\n");

    for (const auto& contig : contigsWithSizes) {
        fprintf(writer_, "##contig=<ID=%s,length=%d>\n", contig.first.c_str(), contig.second);
    }

    fprintf(writer_, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw depth of coverage\">\n");
    fprintf(writer_, "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS mapping quality\">\n");
    fprintf(writer_, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    fprintf(writer_, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
    fprintf(writer_, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n");

    fprintf(writer_, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n");
}

void Vcf::writeRecord(const GenomeRegion& region, int index,
                      bool /* embedded */, bool indelOkArg) {
    if (!writer_) return;

    const PileUp& pu = region.pileUpRegion(index);
    auto bc = pu.baseCall();

    if (!bc.called()) return;

    int pos = region.locus(index) + 1; // VCF is 1-based
    char refBase = region.refBase(index);

    std::string call = bc.callString(indelOkArg);
    if (call.empty() || call == std::string(1, refBase)) return;

    // Determine REF and ALT
    std::string refStr(1, refBase);
    std::string altStr = call;

    // Handle indels
    if (bc.isInsertion()) {
        refStr = std::string(1, refBase);
        altStr = refBase + call;
    } else if (bc.isDeletion()) {
        refStr = refBase + call;
        altStr = std::string(1, refBase);
    }

    int qual = bc.q() * static_cast<int>(pu.count());
    std::string gt = bc.homo ? "1/1" : "0/1";

    fprintf(writer_, "%s\t%d\t.\t%s\t%s\t%d\tPASS\tDP=%lld;MQ=%.1f\tGT:GQ:DP\t%d:%d:%lld\n",
            region.name.c_str(), pos, refStr.c_str(), altStr.c_str(),
            qual, pu.depth(), static_cast<double>(pu.meanMq()),
            bc.homo ? 2 : 1, qual, pu.count());
}

void Vcf::writeFixRecord(const GenomeRegion& region, const GenomeRegion::Fix& fix) {
    if (!writer_) return;

    int pos = std::get<0>(fix) + 1; // 1-based
    const std::string& ref = std::get<1>(fix);
    const std::string& alt = std::get<2>(fix);

    fprintf(writer_, "%s\t%d\t.\t%s\t%s\t.\tPASS\t.\t.\n",
            region.name.c_str(), pos, ref.c_str(), alt.c_str());
}

void Vcf::close() {
    if (writer_) {
        fclose(writer_);
        writer_ = nullptr;
    }
}

} // namespace pilon
