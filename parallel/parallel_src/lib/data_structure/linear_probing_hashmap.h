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
                m_used.resize(m_real_size, NOT_CONTAINED);

                CONTAINED_FLAG = NOT_CONTAINED - 1;
                m_last_request = NOT_CONTAINED;
                m_last_pos     = NOT_CONTAINED;
        }

        void clear() {
                --CONTAINED_FLAG;

                m_last_request = NOT_CONTAINED;
        }

        // find table position or the next free table position if it is not contained
        NodeID contains( NodeID node ) {
                if( m_last_request == node ) return true;
                NodeID hash_value = hash(node);
                for( NodeID i = hash_value; i < m_real_size; ++i) {
                        if( m_internal_map[i].key == node || m_used[i] != CONTAINED_FLAG ) {
                                hash_value = i;
                                break;
                        }
                }

                return m_internal_map[hash_value].key == node && m_used[hash_value] == CONTAINED_FLAG;
        }

        // find table position or the next free table position if it is not contained
        inline NodeID find( NodeID node ) {
                if( m_last_request == node ) return m_last_pos;

                NodeID hash_value = hash(node);
                for( NodeID i = hash_value; i < m_real_size; ++i) {
                        if( m_internal_map[i].key == node || m_used[i] != CONTAINED_FLAG ) {
                                hash_value = i;
                                break;
                        }
                }

                if( m_used[hash_value] != CONTAINED_FLAG ) {
                        m_used[hash_value] = CONTAINED_FLAG;
                        m_internal_map[hash_value].key = node;
                        m_internal_map[hash_value].value = 0;
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
        NodeID CONTAINED_FLAG = NOT_CONTAINED;

        NodeID m_size; 
        NodeID m_real_size; 
        NodeID m_last_pos;
        NodeID m_last_request;
        std::vector< KeyValuePair > m_internal_map;
        std::vector< NodeID > m_used;
};


#endif /* end of include guard: LINEAR_PROBING_HASHMAP_KQ738TKS */
