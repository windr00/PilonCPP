#ifndef PILON_VCF_H
#define PILON_VCF_H

#include <string>
#include <vector>
#include <cstdio>
#include <ctime>
#include "genome.h"

namespace pilon {

class Vcf {
public:
    Vcf(const std::string& filePath,
        const std::vector<std::pair<std::string, int>>& contigsWithSizes = {});

    ~Vcf();

    void writeHeader(const std::vector<std::pair<std::string, int>>& contigsWithSizes);

    void writeRecord(const GenomeRegion& region, int index,
                     bool embedded = false, bool indelOk = true);

    void writeFixRecord(const GenomeRegion& region, const GenomeRegion::Fix& fix);

    void writeDup(const GenomeRegion& region, const Region& dup);

    void close();

    FILE* writer() { return writer_; }

    static void writeChangesFile(const std::vector<GenomeRegion>& regions,
                                 const std::string& changesPath);

    static bool fixContainsN(const GenomeRegion::Fix& fix);

private:
    FILE* writer_;
    std::string filePath_;
};

} // namespace pilon

#endif
