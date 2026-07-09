#pragma once
#include <cstdint>
#include <string>

namespace champhash {

// =============================================================================
// Java-style String.hashCode() — used by Scala for element hashing.
// s[0]*31^(n-1) + s[1]*31^(n-2) + ... + s[n-1], 32-bit wrapping.
// =============================================================================
inline int32_t javaStringHashCode(const std::string& s) {
    int32_t h = 0;
    for (unsigned char c : s) {
        h = 31 * h + static_cast<int32_t>(c);
    }
    return h;
}

// =============================================================================
// Scala MurmurHash3: mix, finalizeHash, productHash
// Matching scala.util.hashing.MurmurHash3 from Scala 2.13.4
// =============================================================================

namespace murmur3 {

inline int32_t mix(int32_t hash, int32_t data) {
    int32_t h = hash;
    int32_t k = data;
    k *= 0xcc9e2d51;
    k = (k << 15) | (static_cast<uint32_t>(k) >> 17);  // rotateLeft(15)
    k *= 0x1b873593;
    h ^= k;
    h = (h << 13) | (static_cast<uint32_t>(h) >> 19);  // rotateLeft(13)
    h = h * 5 + 0xe6546b64;
    return h;
}

inline int32_t avalanche(int32_t hash) {
    int32_t h = hash;
    h ^= static_cast<uint32_t>(h) >> 16;
    h *= 0x85ebca6b;
    h ^= static_cast<uint32_t>(h) >> 13;
    h *= 0xc2b2ae35;
    h ^= static_cast<uint32_t>(h) >> 16;
    return h;
}

inline int32_t finalizeHash(int32_t hash, int length) {
    return avalanche(hash) ^ length;
}

// Seed for product hashing (Scala constant: 0xcafebabe = -889275714)
constexpr int32_t PRODUCT_SEED = static_cast<int32_t>(0xcafebabe);

// productPrefixHash for Tuple3: Java hashCode of "Tuple3"
inline int32_t tuple3PrefixHash() {
    return javaStringHashCode("Tuple3");
}

// Hash for Tuple3(Int, String, String) — exact match for Scala 2.13.4.
// This is what `element.##` returns for the (loc, was, patch) tuple.
inline int32_t productHash3(int a, const std::string& b, const std::string& c) {
    int32_t h = murmur3::mix(murmur3::PRODUCT_SEED, tuple3PrefixHash());
    h = murmur3::mix(h, a);                         // _1: Int, hashCode = value
    h = murmur3::mix(h, javaStringHashCode(b));     // _2: String
    h = murmur3::mix(h, javaStringHashCode(c));     // _3: String
    return murmur3::finalizeHash(h, 3);
}

} // namespace murmur3

// Convenience alias
inline int32_t scalaTuple3Hash(int a, const std::string& b, const std::string& c) {
    return murmur3::productHash3(a, b, c);
}

} // namespace champhash
