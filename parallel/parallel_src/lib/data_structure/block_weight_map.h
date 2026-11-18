/******************************************************************************
 * block_weight_map.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef BLOCK_WEIGHT_MAP_PARHIP_H
#define BLOCK_WEIGHT_MAP_PARHIP_H

#include <vector>

#include "definitions.h"

// Map block -> weight with open addressing (linear probing, power-of-two capacity).
// A block that is not stored has weight 0, so storing 0 is equivalent to erasing;
// the map is only read by key and never iterated.
class block_weight_map {
public:
        block_weight_map() : m_mask(0), m_count(0) {}

        void reserve(size_t n) {
                size_t cap = 16;
                while( cap < 2 * n ) cap <<= 1;
                if( cap > m_keys.size() ) rehash(cap);
        }

        long get( PartitionID key ) const {
                size_t i = slot(key);
                while( m_keys[i] != empty_key() ) {
                        if( m_keys[i] == key ) return m_values[i];
                        i = (i + 1) & m_mask;
                }
                return 0;
        }

        long & ref( PartitionID key ) {
                if( 2 * (m_count + 1) > m_keys.size() ) rehash(m_keys.size() * 2);
                size_t i = slot(key);
                while( m_keys[i] != empty_key() ) {
                        if( m_keys[i] == key ) return m_values[i];
                        i = (i + 1) & m_mask;
                }
                m_keys[i]   = key;
                m_values[i] = 0;
                m_count++;
                return m_values[i];
        }

private:
        static PartitionID empty_key() { return ~(PartitionID)0; }

        size_t slot( PartitionID key ) const {
                return (size_t)((key * 0x9E3779B97F4A7C15ULL) >> 32) & m_mask;
        }

        void rehash( size_t cap ) {
                std::vector< PartitionID > keys(cap, empty_key());
                std::vector< long > values(cap, 0);
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

        std::vector< PartitionID > m_keys;
        std::vector< long > m_values;
        size_t m_mask;
        size_t m_count;
};

#endif
