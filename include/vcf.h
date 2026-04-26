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

#ifndef PILON_VCF_H
#define PILON_VCF_H

#include <string>
#include <vector>
#include <cstdio>
#include "genome.h"

namespace pilon {

class Vcf {
public:
    Vcf(const std::string& filePath,
        const std::vector<std::pair<std::string, int>>& contigsWithSizes = {});

    ~Vcf();

    void writeHeader(const std::vector<std::pair<std::string, int>>& contigsWithSizes);

    void writeRecord(const GenomeRegion& region, int index,
                     bool embedded = false, bool indelOkArg = true);

    void writeFixRecord(const GenomeRegion& region, const GenomeRegion::Fix& fix);

    void close();

private:
    FILE* writer_;
    std::string filePath_;
};

} // namespace pilon

#endif // PILON_VCF_H
