/******************************************************************************
 * ordered_edge_map.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef ORDERED_EDGE_MAP_H
#define ORDERED_EDGE_MAP_H

#include <vector>

#include "definitions.h"
#include "data_structure/hashed_graph.h"

// Map from an undirected edge {source, target} to its accumulated weight, iterated in exactly the
// order in which a hashed_graph (std::unordered_map with hash_hashed_edge) holding the same insertion
// sequence would be iterated. The contraction iterates these maps, and that order is the edge order
// of the coarse graph, so it must not change.
//
// libstdc++'s unordered_map keeps one singly linked list of all nodes. A new key goes to the front of
// its bucket's run if the bucket is not empty, otherwise to the front of the list; a rehash re-inserts
// the nodes in list order with the same rule into the new bucket array. So for a sequence S inserted
// into B buckets, the list consists of the bucket runs in decreasing order of the time each bucket
// was first used, each run in decreasing order of insertion. order() computes this with a counting
// sort, epoch by epoch, using the same rehash policy object as the container.
//
// The emulation relies on libstdc++ internals, so it is used only for the releases on which
// tests/ordered_edge_map_test.cpp was run against the real container (11, 13 and 15, with g++ and
// with clang++). Every other standard library uses the container itself, so the output equals the
// one of the plain hashed_graph on every toolchain.
#if defined(__GLIBCXX__) && defined(_GLIBCXX_RELEASE) && \
    (_GLIBCXX_RELEASE == 11 || _GLIBCXX_RELEASE == 13 || _GLIBCXX_RELEASE == 15)
#define ORDERED_EDGE_MAP_EMULATES_UNORDERED_MAP 1

// a % d for 64-bit a and d by multiplication (Lemire, Kaser and Kurz, "Faster remainder by direct
// computation", 2019), exact for all 64-bit operands with M = 2^128 / d rounded up
inline ULONG fastmod_u64( ULONG a, __uint128_t M, ULONG d ) {
        __uint128_t low    = M * a;
        __uint128_t bottom = ((low & 0xFFFFFFFFFFFFFFFFULL) * d) >> 64;
        __uint128_t top    = (low >> 64) * d;
        return (ULONG)((bottom + top) >> 64);
}

class ordered_edge_map {
public:
        explicit ordered_edge_map( NodeID k ) : m_k(k), m_mask(0) {}

        // the weight of edge {source, target}; a new edge starts with weight 0 and keeps the
        // orientation of its first insertion (as the key of the unordered_map does)
        NodeWeight & weight( NodeID source, NodeID target ) {
                NodeID lo = source < target ? source : target;
                NodeID hi = source < target ? target : source;
                if( 2 * (m_source.size() + 1) > m_slots.size() ) grow();
                size_t i = slot(lo, hi);
                while( m_slots[i] != EMPTY_SLOT() ) {
                        ULONG id = m_slots[i];
                        if( m_lo[id] == lo && m_hi[id] == hi ) return m_weight[id];
                        i = (i + 1) & m_mask;
                }
                ULONG id = m_source.size();
                m_slots[i] = id;
                insert_into_schedule();
                m_source.push_back(source);
                m_target.push_back(target);
                m_lo.push_back(lo);
                m_hi.push_back(hi);
                m_hash.push_back(lo * m_k + hi); // hash_hashed_edge
                m_weight.push_back(0);
                return m_weight.back();
        }

        size_t size() const { return m_source.size(); }

        void prefetch( NodeID source, NodeID target ) const {
                if( m_slots.empty() ) return;
                NodeID lo = source < target ? source : target;
                NodeID hi = source < target ? target : source;
                __builtin_prefetch(&m_slots[slot(lo, hi)]);
        }

        // calls f(source, target, weight) for every edge in unordered_map iteration order
        template < typename F >
        void for_each( F f ) {
                std::vector< ULONG > order;
                compute_order(order);
                for( size_t i = 0; i < order.size(); i++ ) {
                        ULONG id = order[i];
                        f(m_source[id], m_target[id], m_weight[id]);
                }
        }

private:
        static ULONG EMPTY_SLOT() { return ~0ULL; }

        size_t slot( NodeID lo, NodeID hi ) const {
                return (size_t)(((lo * 0x9E3779B97F4A7C15ULL) ^ (hi * 0xC2B2AE3D27D4EB4FULL)) >> 20) & m_mask;
        }

        void grow() {
                size_t cap = m_slots.empty() ? 1024 : 2 * m_slots.size();
                m_slots.assign(cap, EMPTY_SLOT());
                m_mask = cap - 1;
                for( ULONG id = 0; id < m_source.size(); id++ ) {
                        size_t i = slot(m_lo[id], m_hi[id]);
                        while( m_slots[i] != EMPTY_SLOT() ) i = (i + 1) & m_mask;
                        m_slots[i] = id;
                }
        }

        // the container's rehash decision before inserting element number size()
        void insert_into_schedule() {
                if( m_epoch_start.empty() ) {
                        m_bucket_count = 1; // a default constructed unordered_map has one bucket
                }
                std::pair< bool, std::size_t > r = m_policy._M_need_rehash(m_bucket_count, m_source.size(), 1);
                if( r.first ) {
                        m_bucket_count = r.second;
                        m_epoch_start.push_back(m_source.size());
                        m_epoch_buckets.push_back(m_bucket_count);
                } else if( m_epoch_start.empty() ) {
                        m_epoch_start.push_back(0);
                        m_epoch_buckets.push_back(m_bucket_count);
                }
        }

        void compute_order( std::vector< ULONG > & list ) {
                list.clear();
                std::vector< ULONG > seq;
                std::vector< ULONG > bucket;
                std::vector< ULONG > start;
                std::vector< ULONG > touched;
                for( size_t j = 0; j < m_epoch_start.size(); j++ ) {
                        ULONG begin = m_epoch_start[j];
                        ULONG end   = j + 1 < m_epoch_start.size() ? m_epoch_start[j + 1] : m_source.size();
                        size_t B    = m_epoch_buckets[j];
                        // the rehash re-inserts the current list in its order, then the new
                        // elements of this epoch are inserted
                        seq.swap(list);
                        for( ULONG id = begin; id < end; id++ ) seq.push_back(id);

                        bucket.resize(seq.size());
                        start.assign(B, 0);
                        touched.clear();
                        const __uint128_t M = ~(__uint128_t)0 / B + 1;
                        for( size_t i = 0; i < seq.size(); i++ ) {
                                ULONG b   = fastmod_u64(m_hash[seq[i]], M, B);
                                bucket[i] = b;
                                if( start[b]++ == 0 ) touched.push_back(b);
                        }
                        // bucket runs: the most recently first used bucket comes first
                        ULONG offset = 0;
                        for( size_t t = touched.size(); t-- > 0; ) {
                                ULONG b   = touched[t];
                                ULONG c   = start[b];
                                start[b]  = offset;
                                offset   += c;
                        }
                        // within a run: the most recently inserted element comes first
                        list.resize(seq.size());
                        for( size_t i = seq.size(); i-- > 0; ) {
                                list[start[bucket[i]]++] = seq[i];
                        }
                }
        }

        NodeID m_k;
        size_t m_mask;
        std::vector< ULONG > m_slots;
        std::vector< NodeID > m_source, m_target, m_lo, m_hi;
        std::vector< ULONG > m_hash;
        std::vector< NodeWeight > m_weight;

        std::__detail::_Prime_rehash_policy m_policy;
        size_t m_bucket_count;
        std::vector< ULONG > m_epoch_start;
        std::vector< size_t > m_epoch_buckets;
};

#else

// other standard libraries and releases: the container itself
class ordered_edge_map {
public:
        explicit ordered_edge_map( NodeID k ) : m_k(k) {}

        NodeWeight & weight( NodeID source, NodeID target ) {
                hashed_edge he;
                he.k      = m_k;
                he.source = source;
                he.target = target;
                return m_map[he].weight;
        }

        size_t size() const { return m_map.size(); }

        template < typename F >
        void for_each( F f ) {
                for( hashed_graph::iterator it = m_map.begin(); it != m_map.end(); ++it ) {
                        f(it->first.source, it->first.target, it->second.weight);
                }
        }

private:
        NodeID m_k;
        hashed_graph m_map;
};

#endif

#endif
