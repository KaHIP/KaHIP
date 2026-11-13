/******************************************************************************
 * stamped_node_map.h
 *
 * Lookup-only map from NodeID to a small value with the interface of
 * node_hash_map (operator[], find() against end(), erase(), size()). The
 * values live in per-node arrays taken from a small static pool: an entry is
 * present if its stamp equals the map's generation, so creating a map is O(1)
 * and a lookup is a direct index. Pool stores are reused by later maps; when
 * all are taken, a map falls back to node_hash_map. The maps are never
 * iterated, so results are identical to the hash map.
 *****************************************************************************/

#ifndef STAMPED_NODE_MAP_H
#define STAMPED_NODE_MAP_H

#include <algorithm>
#include <vector>

#include "definitions.h"
#include "data_structure/node_hash_map.h"

template <typename V>
class stamped_node_map {
public:
        stamped_node_map() : m_store(nullptr), m_gen(0), m_size(0) {
                std::vector<store*> & p = pool();
                for(size_t i = 0; i < p.size() && !m_store; i++) {
                        if(!p[i]->busy) m_store = p[i];
                }
                if(!m_store && p.size() < MAX_STORES) {
                        p.push_back(new store());
                        m_store = p.back();
                }
                if(m_store) {
                        m_store->busy = true;
                        if(++m_store->gen == 0) {
                                for(slot & e : m_store->slots) e.stamp = 0;
                                m_store->gen = 1;
                        }
                        m_gen = m_store->gen;
                }
        }

        ~stamped_node_map() { if(m_store) m_store->busy = false; }

        stamped_node_map(stamped_node_map && o) noexcept
                : m_store(o.m_store), m_gen(o.m_gen), m_size(o.m_size), m_fallback(std::move(o.m_fallback)) {
                o.m_store = nullptr; o.m_size = 0;
        }

        stamped_node_map & operator=(stamped_node_map && o) noexcept {
                if(this != &o) {
                        if(m_store) m_store->busy = false;
                        m_store = o.m_store; m_gen = o.m_gen; m_size = o.m_size;
                        m_fallback = std::move(o.m_fallback);
                        o.m_store = nullptr; o.m_size = 0;
                }
                return *this;
        }

        stamped_node_map(const stamped_node_map &) = delete;
        stamped_node_map & operator=(const stamped_node_map &) = delete;

        V * find(NodeID key) {
                if(!m_store) return m_fallback.find(key);
                if(key < m_store->slots.size() && m_store->slots[key].stamp == m_gen) return &m_store->slots[key].val;
                return nullptr;
        }

        V * end() { return nullptr; }

        V & operator[](NodeID key) {
                if(!m_store) return m_fallback[key];
                store & s = *m_store;
                if(key >= s.slots.size()) {
                        size_t cap = std::max((size_t)key + 1, 2 * s.slots.size());
                        s.slots.resize(cap);
                }
                slot & e = s.slots[key];
                if(e.stamp != m_gen) {
                        e.stamp = m_gen;
                        e.val   = V();
                        ++m_size;
                }
                return e.val;
        }

        // the value of a key that is present (no stamp check)
        V & present(NodeID key) { return m_store ? m_store->slots[key].val : m_fallback[key]; }

        size_t erase(NodeID key) {
                if(!m_store) return m_fallback.erase(key);
                if(key < m_store->slots.size() && m_store->slots[key].stamp == m_gen) {
                        m_store->slots[key].stamp = 0;
                        --m_size;
                        return 1;
                }
                return 0;
        }

        size_t size() const { return m_store ? m_size : m_fallback.size(); }

private:
        static const size_t MAX_STORES = 4;

        // stamp and value of a key side by side (one cache line per lookup)
        struct slot {
                unsigned stamp = 0;
                V        val   = V();
        };
        struct store {
                std::vector<slot> slots;
                unsigned gen = 0;
                bool busy    = false;
        };
        static std::vector<store*> & pool() { static std::vector<store*> p; return p; }

        store * m_store;
        unsigned m_gen;
        size_t m_size;
        node_hash_map<V> m_fallback;
};

#endif
