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
#include "utils.h"
#include <iostream>
#include <fstream>
#include <sstream>

namespace pilon {

namespace {
    std::string formatTime() {
        time_t now = time(nullptr);
        struct tm* t = localtime(&now);
        char buf[16];
        strftime(buf, sizeof(buf), "%Y%m%d", t);
        return buf;
    }
}

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
    std::string date = formatTime();
    std::string cmdLine;
    for (const auto& a : Pilon::commandArgs) {
        if (!cmdLine.empty()) cmdLine += " ";
        cmdLine += a;
    }

    fprintf(writer_, "##fileformat=VCFv4.1\n");
    fprintf(writer_, "##fileDate=%s\n", date.c_str());
    fprintf(writer_, "##source=\"PilonCpp 1.0.0\"\n");
    fprintf(writer_, "##PILON=\"%s\"\n", cmdLine.c_str());
    fprintf(writer_, "##reference=%s\n", Pilon::genomePath.c_str());

    for (const auto& c : contigsWithSizes) {
        fprintf(writer_, "##contig=<ID=%s,length=%d>\n", c.first.c_str(), c.second);
    }

    fprintf(writer_, "##FILTER=<ID=LowCov,Description=\"Low Coverage of good reads at location\">\n");
    fprintf(writer_, "##FILTER=<ID=Amb,Description=\"Ambiguous evidence in haploid genome\">\n");
    fprintf(writer_, "##FILTER=<ID=Del,Description=\"This base is in a deletion or change event from another record\">\n");
    fprintf(writer_, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Valid read depth; some reads may have been filtered\">\n");
    fprintf(writer_, "##INFO=<ID=TD,Number=1,Type=Integer,Description=\"Total read depth including bad pairs\">\n");
    fprintf(writer_, "##INFO=<ID=PC,Number=1,Type=Integer,Description=\"Physical coverage of valid inserts across locus\">\n");
    fprintf(writer_, "##INFO=<ID=BQ,Number=1,Type=Integer,Description=\"Mean base quality at locus\">\n");
    fprintf(writer_, "##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Mean read mapping quality at locus\">\n");
    fprintf(writer_, "##INFO=<ID=QD,Number=1,Type=Integer,Description=\"Variant confidence/quality by depth\">\n");
    fprintf(writer_, "##INFO=<ID=BC,Number=4,Type=Integer,Description=\"Count of As, Cs, Gs, Ts at locus\">\n");
    if (Pilon::vcfQE)
        fprintf(writer_, "##INFO=<ID=QE,Number=4,Type=Integer,Description=\"Evidence for As, Cs, Gs, Ts weighted by Q & MQ at locus\">\n");
    else
        fprintf(writer_, "##INFO=<ID=QP,Number=4,Type=Integer,Description=\"Percentage of As, Cs, Gs, Ts weighted by Q & MQ at locus\">\n");
    fprintf(writer_, "##INFO=<ID=IC,Number=1,Type=Integer,Description=\"Number of reads with insertion here\">\n");
    fprintf(writer_, "##INFO=<ID=DC,Number=1,Type=Integer,Description=\"Number of reads with deletion here\">\n");
    fprintf(writer_, "##INFO=<ID=XC,Number=1,Type=Integer,Description=\"Number of reads clipped here\">\n");
    fprintf(writer_, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">\n");
    fprintf(writer_, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Fraction of evidence in support of alternate allele(s)\">\n");
    fprintf(writer_, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
    fprintf(writer_, "##INFO=<ID=SVLEN,Number=.,Type=String,Description=\"Difference in length between REF and ALT alleles\">\n");
    fprintf(writer_, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n");
    fprintf(writer_, "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise change from local reassembly\">\n");
    fprintf(writer_, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    fprintf(writer_, "##FORMAT=<ID=AD,Number=.,Type=String,Description=\"Allelic depths for the ref and alt alleles\">\n");
    fprintf(writer_, "##FORMAT=<ID=DP,Number=1,Type=String,Description=\"Approximate read depth\">\n");
    fprintf(writer_, "##ALT=<ID=DUP,Description=\"Possible segmental duplication\">\n");
    fprintf(writer_, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n");
}

void Vcf::writeRecord(const GenomeRegion& region, int index,
                      bool embedded, bool indelOkArg) {
    if (!writer_) return;

    const PileUp& pu = region.pileUpRegion(index);
    auto bc = pu.baseCall();
    int locus = region.locus(index);
    int depth = static_cast<int>(pu.depth());
    int baseDP = static_cast<int>(bc.baseSum);
    int altBaseDP = static_cast<int>(bc.altBaseSum);

    // Scala: indelOk = indelOkArg && index > 0
    bool indelOk = indelOkArg && index > 0;

    std::string rBase, cBase, callType;
    int refDP = 0, altDP = 0;

    if (indelOk && !embedded && bc.isDeletion()) {
        locus -= 1;
        rBase = std::string(1, region.refBase(locus)) + bc.deletion;
        callType = bc.homoIndel ? "1/1" : "0/1";
        int p = pu.delPct();
        cBase = std::string(1, region.refBase(locus));
        refDP = 100 - p;
        altDP = p;
    } else if (indelOk && !embedded && bc.isInsertion()) {
        locus -= 1;
        rBase = std::string(1, region.refBase(locus));
        callType = bc.homoIndel ? "1/1" : "0/1";
        int p = pu.insPct();
        cBase = rBase + bc.insertion;
        refDP = 100 - p;
        altDP = p;
    } else if (bc.homo) {
        char refB = region.refBase(locus);
        if (refB == bc.base || bc.base == 'N') {
            rBase = std::string(1, refB);
            cBase = std::string(1, bc.base);
            callType = "0/0";
            refDP = baseDP;
            altDP = altBaseDP;
        } else {
            rBase = std::string(1, refB);
            cBase = std::string(1, bc.base);
            callType = "1/1";
            refDP = altBaseDP;
            altDP = baseDP;
        }
    } else {
        char refB = region.refBase(locus);
        if (refB == bc.base) {
            rBase = std::string(1, refB);
            cBase = std::string(1, bc.altBase);
            callType = "0/1";
            refDP = baseDP;
            altDP = altBaseDP;
        } else {
            rBase = std::string(1, refB);
            cBase = std::string(1, bc.base);
            callType = "0/1";
            refDP = altBaseDP;
            altDP = baseDP;
        }
    }

    // Filters
    std::string filter = "PASS";
    if (depth < region.minDepth) filter = "LowCov";
    if (!Pilon::diploid && callType == "0/1") {
        if (filter != "PASS") filter += ";Amb";
        else filter = "Amb";
    }
    if (embedded) {
        if (filter != "PASS") filter += ";Del";
        else filter = "Del";
    }

    // Scala: always write records, use "." for ALT when no change
    std::string cBaseVcf = (cBase == rBase || cBase == "N" || cBase == "." || cBase.empty())
                           ? "." : cBase;

    int ac = 0;
    if (callType == "0/1") ac = 1;
    else if (callType == "1/1") ac = 2;

    double af = 0;
    if (refDP + altDP > 0)
        af = static_cast<double>(altDP) / (refDP + altDP);

    char afBuf[16];
    snprintf(afBuf, sizeof(afBuf), "%.2f", af);

    int callDepth = embedded ? static_cast<int>(pu.count()) : static_cast<int>(pu.depth());

    char infoBuf[1024];
    if (Pilon::vcfQE) {
        snprintf(infoBuf, sizeof(infoBuf),
                 "DP=%d;TD=%lld;BQ=%d;MQ=%d;QD=%d;BC=%s;QE=%s;PC=%lld;IC=%lld;DC=%lld;XC=%lld;AC=%d;AF=%s",
                 callDepth,
                 pu.depth() + pu.badPair,
                 pu.meanQual(),
                 pu.meanMq(),
                 bc.q(),
                 pu.baseCount.toString().c_str(),
                 pu.qualSum.toString().c_str(),
                 pu.physCov,
                 pu.insertions,
                 pu.deletions,
                 pu.clips,
                 ac,
                 afBuf);
    } else {
        snprintf(infoBuf, sizeof(infoBuf),
                 "DP=%d;TD=%lld;BQ=%d;MQ=%d;QD=%d;BC=%s;QP=%s;PC=%lld;IC=%lld;DC=%lld;XC=%lld;AC=%d;AF=%s",
                 callDepth,
                 pu.depth() + pu.badPair,
                 pu.meanQual(),
                 pu.meanMq(),
                 bc.q(),
                 pu.baseCount.toString().c_str(),
                 pu.qualSum.toStringPct().c_str(),
                 pu.physCov,
                 pu.insertions,
                 pu.deletions,
                 pu.clips,
                 ac,
                 afBuf);
    }

    int quality = (indelOk && bc.isDeletion()) ? -1 : bc.score;

    fprintf(writer_, "%s\t%d\t.\t%s\t%s\t%s\t%s\t%s\tGT\t%s\n",
            region.name.c_str(), locus + 1,  // VCF 1-based
            rBase.c_str(), cBaseVcf.c_str(),
            quality < 0 ? "." : std::to_string(quality).c_str(),
            filter.c_str(),
            infoBuf,
            callType.c_str());

    // Scala: if (indelOk && bc.indel && !embedded) writeRecord(region, index, bc.isDeletion && bc.homoIndel, false)
    if (indelOk && bc.indel && !embedded) {
        writeRecord(region, index, bc.isDeletion() && bc.homoIndel, false);
    }
}

void Vcf::writeFixRecord(const GenomeRegion& region, const GenomeRegion::Fix& fix) {
    if (!writer_) return;

    int loc = std::get<0>(fix) - 1;
    char rBase = region.refBase(loc);
    std::string ref = std::string(1, rBase) + std::get<1>(fix);
    std::string alt = std::string(1, rBase) + std::get<2>(fix);

    int svlen = static_cast<int>(alt.length()) - static_cast<int>(ref.length());
    int svend = loc + static_cast<int>(ref.length()) - 1;
    std::string svtype = svlen < 0 ? "DEL" : "INS";

    char buf[512];
    snprintf(buf, sizeof(buf),
             "%s\t%d\t.\t%s\t%s\t.\tPASS\tSVTYPE=%s;SVLEN=%d;END=%d%s\tGT\t1/1",
             region.name.c_str(), loc + 1,  // VCF 1-based
             ref.c_str(), alt.c_str(),
             svtype.c_str(), svlen, svend + 1,  // END is 1-based inclusive
             fixContainsN(fix) ? ";IMPRECISE" : "");

    fprintf(writer_, "%s\n", buf);
}

void Vcf::writeDup(const GenomeRegion& region, const Region& dup) {
    if (!writer_) return;

    int loc = dup.start - 1;
    char rBase = region.refBase(loc);

    fprintf(writer_, "%s\t%d\t.\t%c\t<DUP>\t.\tPASS\tSVTYPE=DUP;SVLEN=%d;END=%d;IMPRECISE\tGT\t./.\n",
            region.name.c_str(), loc + 1,
            rBase, dup.size(), dup.stop);
}

void Vcf::writeChangesFile(const std::vector<GenomeRegion>& regions,
                           const std::string& changesPath) {
    std::ofstream changes(changesPath);
    if (!changes.is_open()) return;

    for (const auto& region : regions) {
        const auto& fxs = region.fixes;
        for (const auto& fix : fxs) {
            int loc = std::get<0>(fix);
            const std::string& from = std::get<1>(fix);
            const std::string& to = std::get<2>(fix);
            std::string fromDisplay = from.empty() ? "." : from;
            std::string toDisplay = to.empty() ? "." : to;
            changes << region.name << "\t" << loc << "\t" << fromDisplay << "\t" << toDisplay << "\n";
        }
    }
    changes.close();
}

bool Vcf::fixContainsN(const GenomeRegion::Fix& fix) {
    return std::get<2>(fix).find('N') != std::string::npos;
}

void Vcf::close() {
    if (writer_) {
        fclose(writer_);
        writer_ = nullptr;
    }
}

} // namespace pilon
