/******************************************************************************
 * communication_graph_search_space.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef COMMUNICATION_SEARCH_SPACE_H49ZQ8A4
#define COMMUNICATION_SEARCH_SPACE_H49ZQ8A4

#include <utility>
#include <unordered_map>
#include "data_structure/graph_access.h"
#include "tools/random_functions.h"
#include "partition_config.h"

extern int global_num_nodes;

namespace std {
        template <>
                struct hash< std::pair<NodeID, NodeID> > {
                        public:
                                size_t operator()(const std::pair<NodeID, NodeID> & x) const throw() {
                                        if( x.first < x.second ) 
                                                return x.second*global_num_nodes + x.first;
                                        else 
                                                return x.first*global_num_nodes + x.second;
                                }
                };
}

class communication_graph_search_space {
        public:
                communication_graph_search_space(PartitionConfig & config, NodeID number_of_nodes);
                virtual ~communication_graph_search_space(); 

		void set_graph_ref( graph_access * C);

                bool done() {
                        return m_unsucc_tries >= (NodeID) m_limit || m_have_to_break;
                }; // are we done?

                void commit_status( bool success ) {
                        if(success) m_unsucc_tries = 0;
                        else m_unsucc_tries++;

                        std::pair< NodeID, NodeID > ret_value = m_list_of_pairs[m_last_pointer];
                        if(!success) {
                                m_pair_active.erase(ret_value);
                        } else {
                                std::queue< NodeID > bfsqueue;
                                std::vector< NodeID > touched_nodes;

                                bfsqueue.push(ret_value.first);
                                bfsqueue.push(ret_value.second);

                                m_deepth[ret_value.first]  = 0;
                                m_deepth[ret_value.second]  = 0;
                                int cur_deepth  = 0;

                                touched_nodes.push_back( ret_value.first );
                                touched_nodes.push_back( ret_value.second );

                                while(!bfsqueue.empty()) {
                                        NodeID node = bfsqueue.front();
                                        bfsqueue.pop(); 

                                        if (m_deepth[node] == cur_deepth) {
                                                cur_deepth++;
                                        }

                                        if( cur_deepth > m_search_deepth ) {
                                                break;
                                        }

                                        forall_out_edges((*C), e, node) {
                                                NodeID target = C->getEdgeTarget(e);
                                                if(m_deepth[target] == -1) {
                                                        m_deepth[target] = cur_deepth;
                                                        bfsqueue.push(target);
                                                        touched_nodes.push_back(target);

                                                        std::pair< NodeID, NodeID > p1(ret_value.first, target);
                                                        std::pair< NodeID, NodeID > p2(ret_value.second, target);

                                                        m_pair_active[p1] = true;
                                                        m_pair_active[p2] = true;
                                                }
                                        } endfor
                                }
                                for( unsigned i = 0; i < touched_nodes.size(); i++) {
                                        m_deepth[touched_nodes[i]] = -1;
                                }
                        }
                }

                std::pair< NodeID, NodeID > nextPair() {
                        int starting_pos = m_pointer;
                        m_last_pointer   = m_pointer;

                        std::pair< NodeID, NodeID > ret_value = m_list_of_pairs[m_pointer++];
                        m_pointer = m_pointer == m_limit ? 0 : m_pointer;

			while( (m_pair_active.find(ret_value) ==  m_pair_active.end()) && m_pointer != starting_pos) {
				m_last_pointer = m_pointer;
				ret_value      = m_list_of_pairs[m_pointer++];
				m_pointer      = m_pointer == m_limit ? 0 : m_pointer;
			}
			if( m_pointer == starting_pos ) {
				m_have_to_break = true;
			}
			return ret_value;
                }

        private:
                std::vector< std::pair< NodeID, NodeID > > m_list_of_pairs;                 
                std::unordered_map< std::pair< NodeID, NodeID>, bool > m_pair_active;                 
                std::vector< int > m_deepth;
                int m_limit;
                int m_pointer;
                int m_last_pointer;
                int m_search_deepth;//deepth of the neighborhood
                bool m_have_to_break;
                NodeID m_unsucc_tries;

		PartitionConfig config;
                graph_access * C;
};


#endif /* end of include guard: COMMUNICATION_SEARCH_SPACE_H49ZQ8A4 */
