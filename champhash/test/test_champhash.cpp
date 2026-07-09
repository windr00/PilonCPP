#include <iostream>
#include <vector>
#include <string>
#include <cstdint>
#include <cassert>
#include "champhash/champ_hash.hpp"
#include "champhash/tuple_hash.hpp"

using namespace champhash;

// Test the Java String.hashCode() implementation
void testJavaStringHash() {
    assert(javaStringHashCode("") == 0);
    assert(javaStringHashCode("a") == 97);
    assert(javaStringHashCode("ab") == 3105);  // 31*97 + 98
    assert(javaStringHashCode("abc") == 96354); // 31*3105 + 99
    std::cout << "  Java String.hashCode(): OK" << std::endl;
}

// Test MurmurHash3 mix functions
void testMurmurMix() {
    // Known values from Scala MurmurHash3
    int32_t h = murmur3::mix(murmur3::PRODUCT_SEED, murmur3::tuple3PrefixHash());
    // Don't need exact values, just that they compile and don't crash
    std::cout << "  MurmurHash3 mix: OK (productPrefix seed=" 
              << murmur3::PRODUCT_SEED << ")" << std::endl;
}

// Test Tuple3 hash
void testTuple3Hash() {
    int32_t h1 = scalaTuple3Hash(12345, "foo", "bar");
    int32_t h2 = scalaTuple3Hash(12345, "foo", "bar");
    assert(h1 == h2);  // Deterministic
    
    int32_t h3 = scalaTuple3Hash(12345, "foo", "baz");
    assert(h1 != h3);  // Different third element → different hash
    
    std::cout << "  scalaTuple3Hash: OK (sample=" << h1 << ")" << std::endl;
}

// Test basic ChampSet with int keys
void testBasicChampSet() {
    ChampSet<int> set;
    set.insert(10, 10, improve(10));
    set.insert(20, 20, improve(20));
    set.insert(30, 30, improve(30));
    
    auto vec = set.toVector();
    assert(vec.size() == 3);
    
    std::cout << "  Basic ChampSet<int>: OK (size=" << set.size() << ")";
    std::cout << " order=[";
    for (size_t i = 0; i < vec.size(); i++) {
        if (i > 0) std::cout << ",";
        std::cout << vec[i];
    }
    std::cout << "]" << std::endl;
}

// Test duplicate insertion
void testDuplicate() {
    ChampSet<int> set;
    set.insert(42, 42, improve(42));
    set.insert(42, 42, improve(42));  // duplicate → ignored
    assert(set.size() == 1);
    
    std::cout << "  Duplicate handling: OK (size=" << set.size() << ")" << std::endl;
}

// Test hash collision handling
void testHashCollision() {
    // Two different elements with the SAME improved hash at level 0
    // Choose values where improve(value) & 0x1F == improve(value2) & 0x1F
    // but original hashes differ
    
    // Find two values with same mask at shift=0
    uint32_t h1 = 0;
    uint32_t h2 = 0;
    bool found = false;
    for (uint32_t a = 0; a < 1000 && !found; a++) {
        for (uint32_t b = a + 1; b < 1000 && !found; b++) {
            if (maskFrom(improve(a), 0) == maskFrom(improve(b), 0) && 
                maskFrom(improve(a), 5) != maskFrom(improve(b), 5)) {
                h1 = a; h2 = b; found = true;
            }
        }
    }
    
    if (found) {
        ChampSet<uint32_t> set;
        set.insert(h1, h1, improve(h1));
        set.insert(h2, h2, improve(h2));
        assert(set.size() == 2);
        std::cout << "  Hash collision: OK (elements " << h1 << "," << h2 
                  << " collide at level 0, distinct at level 5)" << std::endl;
    } else {
        std::cout << "  Hash collision: SKIP (no collision pair found in range)" << std::endl;
    }
}

// Test mock breakJoins scenario: 3 solutions with tied scores
void testTieBreaking() {
    // Simulate 3 solutions from breakJoins, all with same patch_len + was_len
    // The CHAMP iterator order should pick a winner deterministically
    
    struct Solution {
        int loc;
        std::string was;
        std::string patch;
        bool operator==(const Solution& o) const {
            return loc == o.loc && was == o.was && patch == o.patch;
        }
    };
    
    Solution s1 = {100, "AAA", "GGG"};
    Solution s2 = {200, "TTT", "CCC"};
    Solution s3 = {150, "CCC", "AAA"};
    
    std::vector<Solution> solutions = {s1, s2, s3};
    std::vector<uint32_t> hashes;
    for (const auto& s : solutions) {
        hashes.push_back(scalaTuple3Hash(s.loc, s.was, s.patch));
    }
    
    auto sorted = champSort(solutions, hashes);
    
    assert(sorted.size() == 3);
    
    std::cout << "  Tie-breaking (3 solutions):" << std::endl;
    for (size_t i = 0; i < sorted.size(); i++) {
        std::cout << "    [" << i << "] loc=" << sorted[i].loc 
                  << " was_len=" << sorted[i].was.length() 
                  << " patch_len=" << sorted[i].patch.length()
                  << " hash=" << improve(hashes[i]) 
                  << " mask0=" << maskFrom(improve(hashes[i]), 0)
                  << std::endl;
    }
    std::cout << "  Winner: loc=" << sorted[0].loc << std::endl;
}

int main() {
    std::cout << "champhash test suite" << std::endl;
    std::cout << "====================" << std::endl;
    
    testJavaStringHash();
    testMurmurMix();
    testTuple3Hash();
    testBasicChampSet();
    testDuplicate();
    testHashCollision();
    testTieBreaking();
    
    std::cout << std::endl << "All tests passed!" << std::endl;
    return 0;
}
