#ifndef PILON_GAPFILLER_H
#define PILON_GAPFILLER_H

#include <string>
#include <vector>
#include <utility>
#include <tuple>
#include "assembler.h"
#include "genome.h"

namespace pilon {

class BamFile;

class GapFiller {
public:
    static constexpr int minExtend = 20;
    static int k;

    GapFiller(GenomeRegion& region);

    std::tuple<int, std::string, std::string> fillGap(const Region& gap);
    std::tuple<int, std::string, std::string> fixBreak(const Region& brk);

    const std::string& tandemRepeat() const { return tandemRepeat_; }
    int breakRadius() const;

    // Static wrappers — return fixes to be applied by fixIssues
    static std::tuple<int, std::string, std::string>
    doFixGap(GenomeRegion& region, const Region& gap);
    static std::tuple<int, std::string, std::string>
    doFixBreak(GenomeRegion& region, const Region& brk);
    static std::vector<std::tuple<int, std::string, std::string>>
    doCloseCircle(GenomeRegion& region, int estimatedLength);
    static void fixNovel(GenomeFile* genome, const std::vector<BamFile*>& bamFiles);

private:
    GenomeRegion& region_;
    std::string tandemRepeat_;
    static const std::tuple<int, std::string, std::string> noSolution;

    // Close circle (matching Scala GapFiller.closeCircle)
    std::vector<std::tuple<int, std::string, std::string>>
    closeCircle(int estimatedLength);

    std::tuple<int, std::string, std::string>
    assembleAcrossBreak(const Region& brk, bool isGap);

    std::string consensusFromLeft(const std::vector<std::string>& seqs) const;
    std::string consensusFromRight(const std::vector<std::string>& seqs) const;

    std::tuple<int, std::vector<std::string>, std::vector<std::string>, int, std::string>
    assembleIntoBreak(const Region& brk, const std::vector<BamRead>& reads);

    std::vector<std::tuple<int, std::string, std::string>>
    breakJoins(int start, const std::vector<std::string>& forwardPaths,
               const std::vector<std::string>& reversePaths, int stop);

    std::tuple<int, std::string, std::string>
    joinBreak(int startArg, const std::string& forward,
              const std::string& reverse, int stopArg);

    static std::string properOverlap(const std::string& left, const std::string& right,
                                     int minOverlap);

    std::tuple<int, std::string, std::string>
    trimPatch(int startArg, const std::string& patchArg, int stopArg);

    bool partialMatchesReference(int start, const std::string& fromLeft,
                                 const std::string& fromRight, int stop,
                                 int loopLength);

    std::vector<BamRead> recruitReads(const Region& brk) const;
    std::vector<BamRead> recruitFrags(const Region& reg) const;
    std::vector<BamRead> recruitJumps(const Region& reg) const;
    std::vector<BamRead> recruitUnpaired(const Region& reg) const;
    std::vector<BamRead> recruitReadsOfType(const Region& reg, const std::string& type) const;
    std::vector<BamRead> recruitReadsFromBams(const Region& reg,
                                              const std::vector<BamFile*>& bams) const;
};

} // namespace pilon

#endif
