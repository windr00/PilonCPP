#pragma once
#include <cstdint>
#include <cstddef>
#include <vector>
#include <memory>
#include <utility>
#include <functional>
#include "hashing.hpp"

namespace champhash {

// CHAMP-based map (key-value store) with iteration order matching
// Scala 2.13.4 scala.collection.immutable.HashMap.

template<typename K, typename V, typename Hasher>
class ChampMap {
    struct Node {
        uint32_t dataMap = 0;
        uint32_t nodeMap = 0;
        std::vector<K> keys;
        std::vector<V> values;
        std::vector<uint32_t> originalHashes;
        std::vector<std::unique_ptr<Node>> subnodes;

        // Collision overflow: elements that hash-collide beyond max trie depth
        std::vector<K> collisionKeys;
        std::vector<V> collisionValues;
        std::vector<uint32_t> collisionHashes;

        int dataIndex(uint32_t bitpos) const { return indexFrom(dataMap, bitpos); }
        int nodeIndex(uint32_t bitpos) const { return indexFrom(nodeMap, bitpos); }
    };

    struct Frame { Node* node; int payloadIdx; int subIdx; };

public:
    ChampMap() : size_(0) {}

    bool insert(const K& key, const V& value) {
        uint32_t orig = hasher_(key), imp = improve(orig);
        bool changed = false, updated = false;
        root_ = insertNode(std::move(root_), key, value, orig, imp, 0, changed, updated);
        if (changed && !updated) size_++;
        return changed && !updated;
    }

    const V* find(const K& key) const {
        return findNode(root_.get(), key, hasher_(key), improve(hasher_(key)), 0);
    }

    V* find(const K& key) {
        return const_cast<V*>(
            static_cast<const ChampMap*>(this)->find(key));
    }

    bool contains(const K& key) const { return find(key) != nullptr; }
    int size() const { return size_; }
    bool empty() const { return size_ == 0; }

    const V& at(const K& key) const { return *find(key); }

    template<typename Visitor>
    void forEach(Visitor&& v) const { forEachNode(root_, v); }

    std::vector<K> keys() const {
        std::vector<K> out; out.reserve(size_);
        forEach([&](const K& k, const V&) { out.push_back(k); });
        return out;
    }

    void clear() { root_.reset(); size_ = 0; }

    bool erase(const K& key) {
        if (!find(key)) return false;
        root_ = eraseNode(std::move(root_), key, hasher_(key), improve(hasher_(key)), 0);
        size_--;
        return true;
    }

    class Iterator {
    public:
        Iterator(Node* n) : node_(n), pi_(0), si_(0) { if (node_) advance(); }
        bool hasNext() const { return ok_; }
        const K& key() const { return node_->keys[pi_]; }
        const V& value() const { return node_->values[pi_]; }
        void next() { pi_++; advance(); }
    private:
        void advance() {
            while (node_) {
                if (pi_ < (int)node_->keys.size()) { ok_ = true; return; }
                if (si_ < (int)node_->subnodes.size()) {
                    stack_.push_back({node_, pi_, si_ + 1});
                    node_ = node_->subnodes[si_].get(); pi_ = 0; si_ = 0;
                    continue;
                }
                if (!stack_.empty()) {
                    auto& f = stack_.back();
                    node_ = f.node; pi_ = f.payloadIdx; si_ = f.subIdx;
                    stack_.pop_back(); continue;
                }
                node_ = nullptr;
            }
            ok_ = false;
        }
        Node* node_; int pi_, si_;
        std::vector<Frame> stack_; bool ok_ = false;
    };

    Iterator iterator() const { return Iterator(root_.get()); }

private:
    static std::unique_ptr<Node> insertNode(std::unique_ptr<Node> node,
                                            const K& key, const V& value,
                                            uint32_t oh, uint32_t ih, int shift,
                                            bool& changed, bool& updated) {
        if (!node) {
            auto n = std::make_unique<Node>();
            if (shift >= MAX_TRIE_DEPTH) {
                n->collisionKeys.push_back(key);
                n->collisionValues.push_back(value);
                n->collisionHashes.push_back(oh);
            } else {
                n->dataMap = bitposFrom(maskFrom(ih, shift));
                n->keys.push_back(key); n->values.push_back(value);
                n->originalHashes.push_back(oh);
            }
            changed = true; return n;
        }

        // At max depth, use collision list
        if (shift >= MAX_TRIE_DEPTH) {
            for (size_t i = 0; i < node->collisionKeys.size(); i++) {
                if (node->collisionHashes[i] == oh && node->collisionKeys[i] == key) {
                    node->collisionValues[i] = value;
                    changed = true; updated = true;
                    return node;
                }
            }
            node->collisionKeys.push_back(key);
            node->collisionValues.push_back(value);
            node->collisionHashes.push_back(oh);
            changed = true;
            return node;
        }

        int mask = maskFrom(ih, shift);
        uint32_t bp = bitposFrom(mask);

        if (node->dataMap & bp) {
            int idx = node->dataIndex(bp);
            if (node->originalHashes[idx] == oh && node->keys[idx] == key) {
                node->values[idx] = value; changed = true; updated = true;
                return node;
            }
            auto sub = std::make_unique<Node>();
            bool d = false, d2 = false;
            sub = insertNode(std::move(sub), node->keys[idx], node->values[idx],
                             node->originalHashes[idx], improve(node->originalHashes[idx]),
                             shift + BIT_PARTITION_SIZE, d, d2);
            sub = insertNode(std::move(sub), key, value, oh, ih,
                             shift + BIT_PARTITION_SIZE, d, d2);
            node->dataMap &= ~bp; node->nodeMap |= bp;
            int si = node->nodeIndex(bp);
            node->keys.erase(node->keys.begin() + idx);
            node->values.erase(node->values.begin() + idx);
            node->originalHashes.erase(node->originalHashes.begin() + idx);
            node->subnodes.insert(node->subnodes.begin() + si, std::move(sub));
            changed = true; return node;
        }
        if (node->nodeMap & bp) {
            int idx = node->nodeIndex(bp);
            bool sc = false, su = false;
            node->subnodes[idx] = insertNode(std::move(node->subnodes[idx]),
                                             key, value, oh, ih,
                                             shift + BIT_PARTITION_SIZE, sc, su);
            if (sc) changed = true;
            if (su) updated = true;
            return node;
        }
        node->dataMap |= bp;
        int di = node->dataIndex(bp);
        node->keys.insert(node->keys.begin() + di, key);
        node->values.insert(node->values.begin() + di, value);
        node->originalHashes.insert(node->originalHashes.begin() + di, oh);
        changed = true; return node;
    }

    static const V* findNode(const Node* node, const K& key,
                              uint32_t oh, uint32_t ih, int shift) {
        if (!node) return nullptr;

        // Check collision list first (at max depth all elements are here)
        for (size_t i = 0; i < node->collisionKeys.size(); i++) {
            if (node->collisionHashes[i] == oh && node->collisionKeys[i] == key)
                return &node->collisionValues[i];
        }

        int mask = maskFrom(ih, shift);
        uint32_t bp = bitposFrom(mask);
        if (node->dataMap & bp) {
            int idx = node->dataIndex(bp);
            if (node->originalHashes[idx] == oh && node->keys[idx] == key)
                return &node->values[idx];
            return nullptr;
        }
        if (node->nodeMap & bp) {
            int idx = node->nodeIndex(bp);
            return findNode(node->subnodes[idx].get(), key, oh, ih,
                           shift + BIT_PARTITION_SIZE);
        }
        return nullptr;
    }

    template<typename Visitor>
    static void forEachNode(const std::unique_ptr<Node>& node, Visitor& v) {
        if (!node) return;
        uint32_t comb = node->dataMap | node->nodeMap;
        int pi = 0, si = 0;
        while (comb) {
            uint32_t bp = comb & -comb; comb ^= bp;
            if (node->dataMap & bp) { v(node->keys[pi], node->values[pi]); pi++; }
            if (node->nodeMap & bp) { forEachNode(node->subnodes[si], v); si++; }
        }
        for (size_t i = 0; i < node->collisionKeys.size(); i++) {
            v(node->collisionKeys[i], node->collisionValues[i]);
        }
    }

    static std::unique_ptr<Node> eraseNode(std::unique_ptr<Node> node,
                                           const K& key, uint32_t oh, uint32_t ih, int shift) {
        if (!node) return nullptr;

        // Check collision list
        for (size_t i = 0; i < node->collisionKeys.size(); i++) {
            if (node->collisionHashes[i] == oh && node->collisionKeys[i] == key) {
                node->collisionKeys.erase(node->collisionKeys.begin() + i);
                node->collisionValues.erase(node->collisionValues.begin() + i);
                node->collisionHashes.erase(node->collisionHashes.begin() + i);
                return node;
            }
        }

        int mask = maskFrom(ih, shift);
        uint32_t bp = bitposFrom(mask);

        if (node->dataMap & bp) {
            int idx = node->dataIndex(bp);
            if (node->originalHashes[idx] == oh && node->keys[idx] == key) {
                node->dataMap &= ~bp;
                node->keys.erase(node->keys.begin() + idx);
                node->values.erase(node->values.begin() + idx);
                node->originalHashes.erase(node->originalHashes.begin() + idx);
                return node;
            }
            return node;
        }
        if (node->nodeMap & bp) {
            int idx = node->nodeIndex(bp);
            node->subnodes[idx] = eraseNode(std::move(node->subnodes[idx]),
                                            key, oh, ih, shift + BIT_PARTITION_SIZE);
            if (node->subnodes[idx] &&
                node->subnodes[idx]->dataMap == 0 &&
                node->subnodes[idx]->nodeMap == 0 &&
                node->subnodes[idx]->collisionKeys.empty()) {
                node->nodeMap &= ~bp;
                node->subnodes.erase(node->subnodes.begin() + idx);
            }
            return node;
        }
        return node;
    }

    Hasher hasher_;
    std::unique_ptr<Node> root_;
    int size_;
};

} // namespace champhash
