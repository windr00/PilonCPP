#pragma once
#include <cstdint>
#include <cstddef>

namespace champhash {

// From Scala 2.13.4: scala.collection.Hashing.improve
inline uint32_t improve(uint32_t hcode) {
    uint32_t h = hcode + ~(hcode << 9);
    h = h ^ (h >> 14);
    h = h + (h << 4);
    h = h ^ (h >> 10);
    return h;
}

// From Scala 2.13.4 ChampCommon
constexpr int BIT_PARTITION_SIZE = 5;
constexpr int BIT_PARTITION_MASK = (1 << BIT_PARTITION_SIZE) - 1;  // 0x1F
constexpr int BRANCHING_FACTOR = 1 << BIT_PARTITION_SIZE;  // 32
constexpr int MAX_DEPTH = 7;  // ceil(32 / 5)
constexpr int MAX_TRIE_DEPTH = BIT_PARTITION_SIZE * (MAX_DEPTH - 1);  // 30, beyond which collision list is used

inline int maskFrom(uint32_t hash, int shift) {
    return (hash >> shift) & BIT_PARTITION_MASK;
}

inline uint32_t bitposFrom(int mask) {
    return 1u << mask;
}

// popcount (bitCount) from <bit> or builtin
inline int popcount(uint32_t x) {
    return __builtin_popcount(x);
}

// indexFrom(bitmap, bitpos): number of set bits in bitmap below bitpos
inline int indexFrom(uint32_t bitmap, uint32_t bitpos) {
    return popcount(bitmap & (bitpos - 1));
}

// indexFrom(bitmap, mask, bitpos): Scala override for full bitmap
inline int indexFrom(uint32_t bitmap, int mask, uint32_t bitpos) {
    if (bitmap == 0xFFFFFFFFu) return mask;
    return indexFrom(bitmap, bitpos);
}

} // namespace champhash
