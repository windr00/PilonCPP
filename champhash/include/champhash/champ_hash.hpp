#pragma once
#include <cstdint>
#include <cstddef>
#include <vector>
#include <algorithm>
#include <memory>
#include <functional>
#include <utility>
#include "hashing.hpp"

namespace champhash {

// CHAMP (Compressed Hash-Array Mapped Prefix-tree) implementation
// matching Scala 2.13.4 scala.collection.immutable.HashSet / BitmapIndexedSetNode.
//
// Key: arbitrary user type T.
// Hash: 32-bit improved hash (computed externally). Internally we store
//       the original (unimproved) hash for collision comparisons.
//
// The trie uses 5-bit partitions (32-way branching). Elements are stored
// in payloads/subnodes in bitpos order.

template<typename T>
class ChampSet {
public:
    ChampSet() : size_(0) {}

    // Insert element with its original (unimproved) hash and improved hash.
    // Returns true if inserted, false if duplicate (already present).
    bool insert(const T& element, uint32_t originalHash, uint32_t improvedHash) {
        bool changed = false;
        root_ = insertNode(std::move(root_), element, originalHash, improvedHash, 0, changed);
        if (changed) size_++;
        return changed;
    }

    // Number of elements
    int size() const { return size_; }

    // Visit all elements in CHAMP iteration order (depth-first, bitpos order).
    template<typename Visitor>
    void forEach(Visitor&& visitor) const {
        forEachNode(root_, std::forward<Visitor>(visitor));
    }

    // Return the first element in iteration order (for tie-breaking).
    const T* first() const {
        return firstNode(root_);
    }

    // Collect elements in CHAMP order into a vector.
    std::vector<T> toVector() const {
        std::vector<T> out;
        out.reserve(size_);
        forEach([&out](const T& e) { out.push_back(e); });
        return out;
    }

private:
    struct Node {
        uint32_t dataMap = 0;
        uint32_t nodeMap = 0;
        std::vector<T> payloads;
        std::vector<uint32_t> originalHashes;
        std::vector<std::unique_ptr<Node>> subnodes;
        // Collision overflow: elements that hash-collide beyond max trie depth
        std::vector<T> collisions;
        std::vector<uint32_t> collisionHashes;

        int dataIndex(uint32_t bitpos) const {
            return indexFrom(dataMap, bitpos);
        }

        int nodeIndex(uint32_t bitpos) const {
            return indexFrom(nodeMap, bitpos);
        }
    };

    static std::unique_ptr<Node> insertNode(std::unique_ptr<Node> node,
                                            const T& element,
                                            uint32_t originalHash,
                                            uint32_t improvedHash,
                                            int shift,
                                            bool& changed) {
        if (!node) {
            auto n = std::make_unique<Node>();
            if (shift >= MAX_TRIE_DEPTH) {
                n->collisions.push_back(element);
                n->collisionHashes.push_back(originalHash);
            } else {
                n->dataMap = bitposFrom(maskFrom(improvedHash, shift));
                n->payloads.push_back(element);
                n->originalHashes.push_back(originalHash);
            }
            changed = true;
            return n;
        }

        // At max depth, use collision list
        if (shift >= MAX_TRIE_DEPTH) {
            for (size_t i = 0; i < node->collisions.size(); i++) {
                if (node->collisionHashes[i] == originalHash && node->collisions[i] == element)
                    return node; // duplicate
            }
            node->collisions.push_back(element);
            node->collisionHashes.push_back(originalHash);
            changed = true;
            return node;
        }

        int mask = maskFrom(improvedHash, shift);
        uint32_t bitpos = bitposFrom(mask);

        // Case 1: existing payload at this bit position -> hash collision or duplicate
        if (node->dataMap & bitpos) {
            int idx = node->dataIndex(bitpos);
            uint32_t existingOrigHash = node->originalHashes[idx];
            const T& existingElement = node->payloads[idx];

            if (existingOrigHash == originalHash && existingElement == element) {
                return node; // duplicate, no change
            }

            // Hash collision: migrate both to subnode at next level
            uint32_t existingImprovedHash = improve(existingOrigHash);
            auto sub = buildSubnode(existingElement, existingOrigHash, existingImprovedHash,
                                      element, originalHash, improvedHash,
                                      shift + BIT_PARTITION_SIZE);

            node->dataMap &= ~bitpos;
            node->nodeMap |= bitpos;

            int subIdx = node->nodeIndex(bitpos);
            node->payloads.erase(node->payloads.begin() + idx);
            node->originalHashes.erase(node->originalHashes.begin() + idx);
            node->subnodes.insert(node->subnodes.begin() + subIdx, std::move(sub));

            changed = true;
            return node;
        }

        // Case 2: existing subnode at this bit position — recurse
        if (node->nodeMap & bitpos) {
            int idx = node->nodeIndex(bitpos);
            bool subChanged = false;
            auto newSub = insertNode(std::move(node->subnodes[idx]),
                                     element, originalHash, improvedHash,
                                     shift + BIT_PARTITION_SIZE, subChanged);
            node->subnodes[idx] = std::move(newSub);
            if (subChanged) changed = true;
            return node;
        }

        // Case 3: empty slot — insert as new payload
        node->dataMap |= bitpos;
        int dataIdx = node->dataIndex(bitpos);
        node->payloads.insert(node->payloads.begin() + dataIdx, element);
        node->originalHashes.insert(node->originalHashes.begin() + dataIdx, originalHash);
        changed = true;
        return node;
    }

    // Build a subnode containing two colliding elements
    static std::unique_ptr<Node> buildSubnode(const T& e1, uint32_t origHash1, uint32_t improvedHash1,
                                                const T& e2, uint32_t origHash2, uint32_t improvedHash2,
                                                int shift) {
        auto sub = std::make_unique<Node>();
        bool dummy = false;
        sub = insertNode(std::move(sub), e1, origHash1, improvedHash1, shift, dummy);
        sub = insertNode(std::move(sub), e2, origHash2, improvedHash2, shift, dummy);
        return sub;
    }

    template<typename Visitor>
    static void forEachNode(const std::unique_ptr<Node>& node, Visitor&& visitor) {
        if (!node) return;

        uint32_t combined = node->dataMap | node->nodeMap;
        int payloadIdx = 0;
        int subnodeIdx = 0;

        while (combined) {
            uint32_t bitpos = combined & -combined;
            combined ^= bitpos;

            if (node->dataMap & bitpos) {
                visitor(node->payloads[payloadIdx++]);
            }
            if (node->nodeMap & bitpos) {
                forEachNode(node->subnodes[subnodeIdx++], std::forward<Visitor>(visitor));
            }
        }
        for (size_t i = 0; i < node->collisions.size(); i++) {
            visitor(node->collisions[i]);
        }
    }

    static const T* firstNode(const std::unique_ptr<Node>& node) {
        if (!node) return nullptr;
        if (!node->payloads.empty()) return &node->payloads[0];
        for (auto& sub : node->subnodes) {
            const T* r = firstNode(sub);
            if (r) return r;
        }
        if (!node->collisions.empty()) return &node->collisions[0];
        return nullptr;
    }

    std::unique_ptr<Node> root_;
    int size_;
};

// Convenience: determine the winner among tied elements using CHAMP iteration order.
template<typename T>
T champWinner(const std::vector<T>& elements, const std::vector<uint32_t>& originalHashes) {
    ChampSet<T> set;
    for (size_t i = 0; i < elements.size(); i++) {
        set.insert(elements[i], originalHashes[i], improve(originalHashes[i]));
    }
    auto vec = set.toVector();
    return vec.empty() ? T{} : vec[0];
}

// Return elements sorted by CHAMP iteration order.
template<typename T>
std::vector<T> champSort(const std::vector<T>& elements, const std::vector<uint32_t>& originalHashes) {
    ChampSet<T> set;
    for (size_t i = 0; i < elements.size(); i++) {
        set.insert(elements[i], originalHashes[i], improve(originalHashes[i]));
    }
    return set.toVector();
}

} // namespace champhash
