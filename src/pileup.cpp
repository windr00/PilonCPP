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

#include "pileup.h"
#include "pilon.h"
#include "bases.h"
#include "utils.h"
#include <algorithm>
#include <unordered_map>
#include <sstream>
#include <array>

namespace pilon {

PileUp::PileUp()
    : mqSum(0), qSum(0), physCov(0), insertSize(0), badPair(0),
      deletions(0), delQual(0), insertions(0), insQual(0), clips(0) {}

long long PileUp::count() const { return baseCount.sum(); }
long long PileUp::depth() const { return baseCount.sum() + deletions; }

int PileUp::baseIndex(char c) const {
    static const std::array<int, 256> table = []() {
        std::array<int, 256> t;
        t.fill(-1);
        t['A'] = 0; t['C'] = 1; t['G'] = 2; t['T'] = 3;
        return t;
    }();
    return table[static_cast<unsigned char>(c)];
}

char PileUp::indexBase(int i) const {
    return "ACGT"[i];
}

int PileUp::weightedMq() const { return Utils::roundDiv(qualSum.sum(), qSum); }
int PileUp::weightedQual() const { return Utils::roundDiv(qualSum.sum(), mqSum); }
int PileUp::meanQual() const { return Utils::roundDiv(qualSum.sum(), Utils::roundDiv(mqSum * count(), depth())); }
int PileUp::meanMq() const { return Utils::roundDiv(mqSum - depth(), depth()); }

void PileUp::add(char base, int qual, int mq) {
    int bi = baseIndex(base);
    if (bi >= 0 && qual >= Pilon::minQual) {
        int mq1 = mq + 1;
        baseCount.add(bi);
        qualSum.add(bi, qual * mq1);
        mqSum += mq1;
        qSum += qual;
    }
}

void PileUp::remove(char base, int qual, int mq) {
    int bi = baseIndex(base);
    if (bi >= 0) {
        int mq1 = mq + 1;
        baseCount.remove(bi);
        qualSum.remove(bi, qual * mq1);
        mqSum -= mq1;
        qSum -= qual;
    }
}

void PileUp::addInsertion(const std::vector<uint8_t>& insertion, int qual, int mq) {
    int mq1 = mq + 1;
    insQual += mq1;
    qSum += qual;
    insertionList.push_back(insertion);
    insertions++;
}

void PileUp::addDeletion(const std::vector<uint8_t>& deletion, int qual, int mq) {
    int mq1 = mq + 1;
    mqSum += mq1;
    delQual += mq1;
    qSum += qual;
    deletionList.push_back(deletion);
    deletions++;
}

long long PileUp::totalQSum() const { return qualSum.sum() + insQual + delQual; }

int PileUp::insPct() const {
    return std::max(Utils::pct(insQual, mqSum), Utils::pct(insertions, count()));
}

int PileUp::delPct() const {
    return std::max(Utils::pct(delQual, mqSum), Utils::pct(deletions, count() + deletions));
}

int PileUp::qualPct() const {
    auto [base, max, sum] = qualSum.maxBase();
    return Utils::pct(max, sum);
}

int PileUp::clipPct() const { return Utils::pct(clips, depth()); }

// Helper functions for BaseCall
static std::string doIndelCall(const std::vector<std::vector<uint8_t>>& indelList, int pct, long long depth) {
    if (depth < Pilon::minMinDepth || pct < 33 || indelList.empty()) return "";
    std::unordered_map<std::string, int> map;
    for (const auto& indel : indelList) {
        std::string indelStr(indel.begin(), indel.end());
        map[indelStr]++;
    }
    auto winner = std::max_element(map.begin(), map.end(),
        [](const auto& a, const auto& b) { return a.second < b.second; });
    if (winner->second < static_cast<int>(indelList.size()) / 2) return "";
    const std::string& winStr = winner->first;
    if (pct >= 50 - static_cast<int>(winStr.size())) return winStr;
    return "";
}

static std::pair<std::string, bool> doHetIndelCall(const std::vector<std::vector<uint8_t>>& indelList, int pct, long long depth) {
    if (depth < Pilon::minMinDepth || pct < 5 || indelList.empty()) return {"", true};

    std::unordered_map<std::string, int> map;
    for (const auto& indel : indelList) {
        std::string indelStr(indel.begin(), indel.end());
        map[indelStr]++;
    }
    auto winner = std::max_element(map.begin(), map.end(),
        [](const auto& a, const auto& b) { return a.second < b.second; });
    if (winner->second < 2 || winner->second <= static_cast<int>(indelList.size()) / 2) return {"", true};
    const std::string& winStr = winner->first;
    if (winStr.find('N') != std::string::npos) return {"", true};

    if (Pilon::oldIndel) {
        if (pct >= 33 && pct >= 50 - static_cast<int>(winStr.size()))
            return {winStr, true};
        return {"", true};
    } else {
        int middle = std::max(45 - static_cast<int>(winStr.size()), 10);
        int low = middle / 2;
        int high = middle + middle - low;
        if (pct > high) return {winStr, true};
        else if (pct >= low) return {winStr, false};
        else return {"", true};
    }
}

// BaseCall implementation
PileUp::BaseCall::BaseCall(const PileUp& pu)
    : n(pu.count()) {
    auto order = pu.qSum > 0 ? pu.qualSum.order() : pu.baseCount.order();
    baseIndex = order[0];
    altBaseIndex = order[1];
    base = n > 0 ? pu.indexBase(baseIndex) : 'N';
    baseSum = pu.qualSum.sums[baseIndex];
    altBase = pu.indexBase(altBaseIndex);
    altBaseSum = pu.qualSum.sums[altBaseIndex];

    long long total = pu.qualSum.sum();
    long long homoScore = baseSum - (total - baseSum);
    long long halfTotal = total / 2;
    long long heteroScore = total - std::abs(halfTotal - baseSum) - std::abs(halfTotal - altBaseSum);
    homo = homoScore >= heteroScore;
    score = pu.mqSum > 0 ? std::abs(homoScore - heteroScore) * n / pu.mqSum : 0;

    // Scala insertCall: insertions > 2 && insertions > deletions
    std::string insert_call;
    bool homoIns = true;
    if (pu.insertions > 2 && pu.insertions > pu.deletions) {
        auto [ins, homo] = doHetIndelCall(pu.insertionList, pu.insPct(), pu.depth());
        insert_call = ins;
        homoIns = homo;
    }

    if (!insert_call.empty()) {
        insertion = insert_call;
        deletion = "";
        indel = true;
        homoIndel = homoIns;
    } else {
        // Scala deletionCall: deletions > 2 && deletions > insertions
        if (pu.deletions > 2 && pu.deletions > pu.insertions) {
            auto [del, homoDel] = doHetIndelCall(pu.deletionList, pu.delPct(), pu.depth());
            if (!del.empty()) {
                insertion = "";
                deletion = del;
                indel = true;
                homoIndel = homoDel;
            } else {
                insertion = "";
                deletion = "";
                indel = false;
                homoIndel = true;
            }
        } else {
            insertion = "";
            deletion = "";
            indel = false;
            homoIndel = true;
        }
    }
}

bool PileUp::BaseCall::isInsertion() const { return !insertion.empty(); }
bool PileUp::BaseCall::isDeletion() const { return !deletion.empty(); }
bool PileUp::BaseCall::called() const { return base != 'N' || indel; }
int PileUp::BaseCall::q() const { return n > 0 ? score / n : 0; }
bool PileUp::BaseCall::highConfidence() const { return q() >= 10; }

std::string PileUp::BaseCall::callString(bool indelOk) const {
    if (indelOk && isInsertion()) return insertion;
    if (indelOk && isDeletion()) return deletion;
    return std::string(1, base);
}

bool PileUp::BaseCall::baseMatch(char refBase) const { return refBase == base; }

char PileUp::BaseCall::iupacBase() const {
    return homo ? base : Bases::toIUPAC(base, altBase);
}

std::string PileUp::BaseCall::insertCall() const {
    return insertion;
}

std::string PileUp::BaseCall::deletionCall() const {
    return deletion;
}

std::string PileUp::BaseCall::indelCall(const std::vector<std::vector<uint8_t>>& indelList, int pct) const {
    return doIndelCall(indelList, pct, n);
}

std::pair<std::string, bool> PileUp::BaseCall::hetIndelCall(
    const std::vector<std::vector<uint8_t>>& indelList, int pct) const {
    return doHetIndelCall(indelList, pct, n);
}

std::string PileUp::BaseCall::toString() const {
    std::ostringstream oss;
    if (isInsertion()) oss << "bc=i" << insertion;
    else if (isDeletion()) oss << "bc=d" << deletion;
    else {
        oss << "bc=" << base;
        if (!homo) oss << "/" << altBase;
        oss << ",cq=" << q();
    }
    return oss.str();
}

PileUp::BaseCall PileUp::baseCall() const { return BaseCall(*this); }

std::string PileUp::toString() const {
    std::ostringstream oss;
    oss << "<PileUp " << baseCall().toString() << ",b=" << baseCount.toString()
        << "/" << qualSum.toStringPct() << ",c=" << depth() << "/" << (depth() + badPair)
        << ",i=" << insertions << ",d=" << deletions << ",q=" << weightedQual()
        << ",mq=" << weightedMq() << ",p=" << physCov << ",s=" << insertSize
        << ",x=" << clips << ">";
    return oss.str();
}

} // namespace pilon
