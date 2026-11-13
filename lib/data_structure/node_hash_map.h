/******************************************************************************
 * node_hash_map.h
 *
 * Open-addressing (linear probing) hash map from NodeID to a small value.
 * Replaces std::unordered_map where the map is only used for lookups, so the
 * iteration order of std::unordered_map never matters. Provides the subset of
 * the std::unordered_map interface used by KaHIP: operator[] (inserts a
 * value-initialized entry when absent), find() compared against end()
 * (pointers, nullptr when absent), erase(), size().
 *****************************************************************************/

#ifndef NODE_HASH_MAP_H
#define NODE_HASH_MAP_H

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>

#include "definitions.h"

template <typename V>
class node_hash_map {
public:
        node_hash_map() : m_size(0), m_mask(0) {}

        V * find(NodeID key) {
                if(m_size == 0) return nullptr;
                size_t i = slot(key);
                for(;;) {
                        NodeID k = m_keys[i];
                        if(k == key) return &m_vals[i];
                        if(k == EMPTY) return nullptr;
                        i = (i + 1) & m_mask;
                }
        }

        V * end() { return nullptr; }

        V & operator[](NodeID key) {
                if(2 * (m_size + 1) > m_keys.size()) grow();
                size_t i = slot(key);
                for(;;) {
                        NodeID k = m_keys[i];
                        if(k == key) return m_vals[i];
                        if(k == EMPTY) {
                                m_keys[i] = key;
                                m_vals[i] = V();
                                ++m_size;
                                return m_vals[i];
                        }
                        i = (i + 1) & m_mask;
                }
        }

        size_t erase(NodeID key) {
                if(m_size == 0) return 0;
                size_t i = slot(key);
                for(;;) {
                        NodeID k = m_keys[i];
                        if(k == EMPTY) return 0;
                        if(k == key) break;
                        i = (i + 1) & m_mask;
                }
                // backward-shift deletion
                size_t j = i;
                for(;;) {
                        j = (j + 1) & m_mask;
                        NodeID k = m_keys[j];
                        if(k == EMPTY) break;
                        size_t home = slot(k);
                        // move k to i if its home is not cyclically in (i, j]
                        bool in_range = (i <= j) ? (i < home && home <= j) : (i < home || home <= j);
                        if(!in_range) {
                                m_keys[i] = k;
                                m_vals[i] = m_vals[j];
                                i = j;
                        }
                }
                m_keys[i] = EMPTY;
                --m_size;
                return 1;
        }

        size_t size() const { return m_size; }

private:
        static const NodeID EMPTY = std::numeric_limits<NodeID>::max();

        size_t slot(NodeID key) const {
                return (size_t)(((uint64_t)key * 0x9E3779B97F4A7C15ULL) >> 32) & m_mask;
        }

        void grow() {
                size_t cap = m_keys.empty() ? 32 : 2 * m_keys.size();
                std::vector<NodeID> old_keys(cap, EMPTY);
                std::vector<V> old_vals(cap);
                old_keys.swap(m_keys);
                old_vals.swap(m_vals);
                m_mask = cap - 1;
                for(size_t i = 0; i < old_keys.size(); i++) {
                        if(old_keys[i] != EMPTY) {
                                size_t j = slot(old_keys[i]);
                                while(m_keys[j] != EMPTY) j = (j + 1) & m_mask;
                                m_keys[j] = old_keys[i];
                                m_vals[j] = old_vals[i];
                        }
                }
        }

        std::vector<NodeID> m_keys;
        std::vector<V> m_vals;
        size_t m_size;
        size_t m_mask;
};

// Out-of-class definition: grow() binds EMPTY to a const reference (odr-use), which needs it in C++11.
template <typename V>
const NodeID node_hash_map<V>::EMPTY;

#endif
