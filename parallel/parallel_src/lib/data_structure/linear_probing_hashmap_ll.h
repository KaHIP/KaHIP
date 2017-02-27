/******************************************************************************
 * linear_probing_hashmap_ll.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Graph Partitioning 
 ******************************************************************************
 * Copyright (C) 2017 Christian Schulz 
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#ifndef LINEAR_PROBING_HASHMAP_LL_KQ738TKS
#define LINEAR_PROBING_HASHMAP_LL_KQ738TKS

#include <stack>

const ULONG NOT_CONTAINED_LL = std::numeric_limits<ULONG>::max();

struct KeyValuePair {
        ULONG key;
        ULONG value; // may be gain or a node id

        KeyValuePair() {
                key = NOT_CONTAINED_LL;
                value = 0;
        }
};

// no resize, no deletes
class linear_probing_hashmap_ll {
public:
        linear_probing_hashmap_ll() {};

        virtual ~linear_probing_hashmap_ll() {};
        void init( ULONG max_size, ULONG factor ) {
                m_size      = max_size*factor;
                m_real_size = m_size + max_size*1.1; // allow some slack in the end so that we do not need to handle special cases

                m_internal_map.resize(m_real_size);
                m_last_request = NOT_CONTAINED_LL;
                m_last_pos     = NOT_CONTAINED_LL;
        }

        void clear() {
                while( m_contained_key_positions.size() > 0 ) {
                        ULONG cur_key_position = m_contained_key_positions.top();
                        m_contained_key_positions.pop();

                        m_internal_map[cur_key_position].key   = NOT_CONTAINED_LL;
                        m_internal_map[cur_key_position].value = 0;
                }
                m_last_request = NOT_CONTAINED_LL;
        }

        // find table position or the next free table position if it is not contained
        ULONG contains( ULONG node ) {
                if( m_last_request == node ) return true;
                ULONG hash_value = hash(node);
                for( ULONG i = hash_value; i < m_real_size; i++) {
                        if( m_internal_map[i].key == node || m_internal_map[i].key == NOT_CONTAINED_LL) {
                                hash_value = i;
                                break;
                        }
                }

                return m_internal_map[hash_value].key == node;
        }

        // find table position or the next free table position if it is not contained
        ULONG find( ULONG node ) {
                if( m_last_request == node ) return m_last_pos;

                ULONG hash_value = hash(node);
                for( ULONG i = hash_value; i < m_real_size; i++) {
                        if( m_internal_map[i].key == node || m_internal_map[i].key == NOT_CONTAINED_LL) {
                                hash_value = i;
                                break;
                        }
                }

                if( m_internal_map[hash_value].key == NOT_CONTAINED_LL ) {
                        m_internal_map[hash_value].key = node;
                        m_internal_map[hash_value].value = 0;
                        m_contained_key_positions.push(hash_value);
                }

                m_last_request = node;
                m_last_pos     = hash_value;
                return hash_value;
        }

        ULONG& operator[](ULONG node) {
                ULONG table_pos = find(node);
                return m_internal_map[table_pos].value;
        };

        ULONG hash(ULONG key) {
                return key % m_size;
        } 

        std::vector< KeyValuePair > * internal_access() {return & m_internal_map; }
private:
        ULONG m_size; 
        ULONG m_real_size; 
        ULONG m_last_pos;
        ULONG m_last_request;
        std::vector< KeyValuePair > m_internal_map;
        std::stack< ULONG > m_contained_key_positions;
};


#endif /* end of include guard: LINEAR_PROBING_HASHMAP_KQ738TKS */
