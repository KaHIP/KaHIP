/******************************************************************************
 * linear_probing_hashmap.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef LINEAR_PROBING_HASHMAP_KQ738TKS
#define LINEAR_PROBING_HASHMAP_KQ738TKS

#include <vector>

const NodeID NOT_CONTAINED = std::numeric_limits<NodeID>::max();

struct KeyValuePair {
        NodeID key;
        NodeID value; // may be gain or a node id

        KeyValuePair() {
                key = NOT_CONTAINED;
                value = 0;
        }
};

// no resize, no deletes
// Open addressing with linear probing in a power-of-two table that is at most half full
// (multiplicative hash); clear() resets only the slots used since the last clear.
class linear_probing_hashmap {
public:
        linear_probing_hashmap() {};
        virtual ~linear_probing_hashmap() {};

        void init( ULONG max_size, ULONG factor ) {
                size_t cap = 16;
                while( cap < 2 * (size_t)max_size ) cap <<= 1;
                m_internal_map.assign(cap, KeyValuePair());
                m_mask         = cap - 1;
                m_contained_key_positions.clear();
                m_contained_key_positions.reserve(max_size);
                m_last_request = NOT_CONTAINED;
                m_last_pos     = NOT_CONTAINED;
        }

        // use only the first part of the table that keeps a node of this degree at most half full;
        // only valid while the map is empty (between two clear() calls the keys stay in place)
        void use_capacity_for( NodeID degree ) {
                size_t cap = 16;
                while( cap < 2 * (size_t)degree ) cap <<= 1;
                if( cap > m_internal_map.size() ) cap = m_internal_map.size();
                m_mask = cap - 1;
        }

        void clear() {
                for( size_t i = 0; i < m_contained_key_positions.size(); i++) {
                        KeyValuePair & kv = m_internal_map[m_contained_key_positions[i]];
                        kv.key   = NOT_CONTAINED;
                        kv.value = 0;
                }
                m_contained_key_positions.clear();
                m_last_request = NOT_CONTAINED;
        }

        // find table position or the next free table position if it is not contained
        NodeID contains( NodeID node ) {
                if( m_last_request == node ) return true;
                NodeID pos = hash(node);
                while( m_internal_map[pos].key != node && m_internal_map[pos].key != NOT_CONTAINED ) {
                        pos = (pos + 1) & m_mask;
                }
                return m_internal_map[pos].key == node;
        }

        // find table position or the next free table position if it is not contained
        NodeID find( NodeID node ) {
                if( m_last_request == node ) return m_last_pos;

                NodeID pos = hash(node);
                while( m_internal_map[pos].key != node && m_internal_map[pos].key != NOT_CONTAINED ) {
                        pos = (pos + 1) & m_mask;
                }

                if( m_internal_map[pos].key == NOT_CONTAINED ) {
                        m_internal_map[pos].key   = node;
                        m_internal_map[pos].value = 0;
                        m_contained_key_positions.push_back(pos);
                }

                m_last_request = node;
                m_last_pos     = pos;
                return pos;
        }

        NodeID& operator[](NodeID node) {
                NodeID table_pos = find(node);
                return m_internal_map[table_pos].value;
        };

        NodeID hash(NodeID key) {
                return (NodeID)((key * 0x9E3779B97F4A7C15ULL) >> 32) & m_mask;
        } 

        std::vector< KeyValuePair > * internal_access() {return & m_internal_map; }

private:
        NodeID m_mask;
        NodeID m_last_pos;
        NodeID m_last_request;
        std::vector< KeyValuePair > m_internal_map;
        std::vector< NodeID > m_contained_key_positions;
};


#endif /* end of include guard: LINEAR_PROBING_HASHMAP_KQ738TKS */
