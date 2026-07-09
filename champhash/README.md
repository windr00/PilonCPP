# champhash

C++17 header-only library implementing Scala 2.13.4 `scala.collection.immutable.HashSet`
iteration order, based on the **CHAMP** (Compressed Hash-Array Mapped Prefix-tree) algorithm.

## Purpose

Scala's immutable `HashSet` traverses elements in a deterministic order governed by
the CHAMP trie structure (5-bit partitions, 32-way branching). This library replicates
that exact ordering in C++, enabling byte-level compatibility with Scala's internal
HashSet iteration order.

## Files

| File | Description |
|------|-------------|
| `include/champhash/hashing.hpp` | CHAMP constants (`BitPartitionSize=5`), `improve()`, `maskFrom()`, `bitposFrom()`, `indexFrom()` |
| `include/champhash/champ_hash.hpp` | `ChampSet<T>` — CHAMP trie container |
| `include/champhash/tuple_hash.hpp` | Scala `MurmurHash3.productHash` for `Tuple3(Int,String,String)` |

## Building

```bash
mkdir build && cd build
cmake .. -DBUILD_TESTS=ON
make -j$(nproc)
./test/test_champhash
```

## API

### `ChampSet<T>`

```cpp
ChampSet<T> set;

// Insert with original hash and improved hash
bool inserted = set.insert(element, originalHash, improve(originalHash));

// Iterate in CHAMP order (matching Scala HashSet)
set.forEach([&](const T& e) { ... });

// Get first element (winner in tie-breaking)
const T* winner = set.first();

// Collect all elements in CHAMP order
std::vector<T> vec = set.toVector();
```

### Convenience functions

```cpp
// Sort by CHAMP iteration order
auto sorted = champSort(elements, originalHashes);

// Pick the winner (first in CHAMP order)
auto winner = champWinner(elements, originalHashes);
```

## Design notes

- **Not a general-purpose hash set**: optimized for < 12 elements per node,
  matching break-fix scenarios. Uses `std::vector` for low overhead.
- **Iterates in CHAMP order**: dataMap bits (LSB → MSB) first, then subnodes
  recursively in nodeMap bit order.
- **Immutable-style insertion**: `insertNode` returns a new root pointer, though
  mutates in place for efficiency (no structural sharing needed).
- **Hash collision handling**: when two elements collide at the same 5-bit
  partition, a subnode is created at the next shift level (`shift + 5`).

## Reference

- Scala 2.13.4 `scala/collection/immutable/HashSet.scala`
- Scala 2.13.4 `scala/collection/ChampCommon.scala`
- M. J. Steindorfer and J. J. Vinju, "Optimizing Hash-Array Mapped Tries for
  Fast and Lean Immutable JVM Collections" (OOPSLA 2015)
