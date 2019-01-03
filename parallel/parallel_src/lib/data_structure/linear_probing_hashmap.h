/******************************************************************************
 * linear_probing_hashmap.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef LINEAR_PROBING_HASHMAP_KQ738TKS
#define LINEAR_PROBING_HASHMAP_KQ738TKS

#include <stack>

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
class linear_probing_hashmap {
public:
        linear_probing_hashmap() {};

        virtual ~linear_probing_hashmap() {};
        void init( ULONG max_size, ULONG factor ) {
                m_size      = max_size*factor;
                m_real_size = m_size + max_size*1.1; // allow some slack in the end so that we do not need to handle special cases

                m_internal_map.resize(m_real_size);
                m_last_request = NOT_CONTAINED;
                m_last_pos     = NOT_CONTAINED;
        }

        void clear() {
                while( m_contained_key_positions.size() > 0 ) {
                        NodeID cur_key_position = m_contained_key_positions.top();
                        m_contained_key_positions.pop();

                        m_internal_map[cur_key_position].key   = NOT_CONTAINED;
                        m_internal_map[cur_key_position].value = 0;
                }
                m_last_request = NOT_CONTAINED;
        }

        // find table position or the next free table position if it is not contained
        NodeID contains( NodeID node ) {
                if( m_last_request == node ) return true;
                NodeID hash_value = hash(node);
                for( NodeID i = hash_value; i < m_real_size; i++) {
                        if( m_internal_map[i].key == node || m_internal_map[i].key == NOT_CONTAINED) {
                                hash_value = i;
                                break;
                        }
                }

                return m_internal_map[hash_value].key == node;
        }

        // find table position or the next free table position if it is not contained
        NodeID find( NodeID node ) {
                if( m_last_request == node ) return m_last_pos;

                NodeID hash_value = hash(node);
                for( NodeID i = hash_value; i < m_real_size; i++) {
                        if( m_internal_map[i].key == node || m_internal_map[i].key == NOT_CONTAINED) {
                                hash_value = i;
                                break;
                        }
                }

                if( m_internal_map[hash_value].key == NOT_CONTAINED ) {
                        m_internal_map[hash_value].key = node;
                        m_internal_map[hash_value].value = 0;
                        m_contained_key_positions.push(hash_value);
                }

                m_last_request = node;
                m_last_pos     = hash_value;
                return hash_value;
        }

        NodeID& operator[](NodeID node) {
                NodeID table_pos = find(node);
                return m_internal_map[table_pos].value;
        };

        NodeID hash(NodeID key) {
                return key % m_size;
        } 

        std::vector< KeyValuePair > * internal_access() {return & m_internal_map; }
private:
        NodeID m_size; 
        NodeID m_real_size; 
        NodeID m_last_pos;
        NodeID m_last_request;
        std::vector< KeyValuePair > m_internal_map;
        std::stack< NodeID > m_contained_key_positions;
};


#endif /* end of include guard: LINEAR_PROBING_HASHMAP_KQ738TKS */
