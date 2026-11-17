/******************************************************************************
 * open_hash_map.h
 * *
 * Map from NodeID to a value with open addressing (linear probing, power-of-two
 * capacity, at most half full, multiplicative hash). The largest NodeID marks an
 * empty slot and cannot be a key. No deletes. Iteration order is deterministic
 * but differs from std::unordered_map, so only use it where the order of
 * for_each does not influence any result.
 *****************************************************************************/

#ifndef OPEN_HASH_MAP_PARHIP_H
#define OPEN_HASH_MAP_PARHIP_H

#include <limits>
#include <vector>

#include "definitions.h"

template < typename V >
class open_hash_map {
public:
        open_hash_map() : m_mask(0), m_count(0) {}

        void reserve( size_t n ) {
                size_t cap = 16;
                while( cap < 2 * n ) cap <<= 1;
                if( cap > m_keys.size() ) rehash(cap);
        }

        // value of key, inserted with V() if it is not contained
        V & ref( NodeID key ) {
                if( 2 * (m_count + 1) > m_keys.size() ) rehash(m_keys.empty() ? 16 : 2 * m_keys.size());
                size_t i = slot(key);
                while( m_keys[i] != empty_key() ) {
                        if( m_keys[i] == key ) return m_values[i];
                        i = (i + 1) & m_mask;
                }
                m_keys[i]   = key;
                m_values[i] = V();
                m_count++;
                return m_values[i];
        }

        // pointer to the value of key, or NULL if it is not contained
        V * find( NodeID key ) {
                if( m_keys.empty() ) return NULL;
                size_t i = slot(key);
                while( m_keys[i] != empty_key() ) {
                        if( m_keys[i] == key ) return &m_values[i];
                        i = (i + 1) & m_mask;
                }
                return NULL;
        }

        size_t size() const { return m_count; }

        template < typename F >
        void for_each( F f ) const {
                for( size_t i = 0; i < m_keys.size(); i++) {
                        if( m_keys[i] != empty_key() ) f(m_keys[i], m_values[i]);
                }
        }

private:
        static NodeID empty_key() { return std::numeric_limits< NodeID >::max(); }

        size_t slot( NodeID key ) const {
                return (size_t)((key * 0x9E3779B97F4A7C15ULL) >> 32) & m_mask;
        }

        void rehash( size_t cap ) {
                std::vector< NodeID > keys(cap, empty_key());
                std::vector< V > values(cap);
                keys.swap(m_keys);
                values.swap(m_values);
                m_mask = cap - 1;
                for( size_t j = 0; j < keys.size(); j++) {
                        if( keys[j] == empty_key() ) continue;
                        size_t i = slot(keys[j]);
                        while( m_keys[i] != empty_key() ) i = (i + 1) & m_mask;
                        m_keys[i]   = keys[j];
                        m_values[i] = values[j];
                }
        }

        std::vector< NodeID > m_keys;
        std::vector< V > m_values;
        size_t m_mask;
        size_t m_count;
};

#endif
