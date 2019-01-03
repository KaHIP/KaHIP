/******************************************************************************
 * full_search_space_pruned.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef FULL_SEARCH_SPACE_PRUNED_H49ZQ8A4
#define FULL_SEARCH_SPACE_PRUNED_H49ZQ8A4

#include <utility>
#include "data_structure/graph_access.h"
#include "partition_config.h"

class full_search_space_pruned {
public:
        full_search_space_pruned(PartitionConfig & config, NodeID number_of_nodes);
        virtual ~full_search_space_pruned();

        void set_graph_ref( graph_access * C ) {}

        bool done() {
             return m_unsucc_tries > m_ub && m_internal_k+1 == ceil(m_number_of_nodes/(double) config.search_space_s);
        }; // are we done?

        void commit_status( bool success ) {
                if(success) m_unsucc_tries = 0;
                else m_unsucc_tries++;
	}

        std::pair< NodeID, NodeID > nextPair() {
                if( m_unsucc_tries > m_ub && m_internal_k+1 != ceil(m_number_of_nodes/(double) config.search_space_s)) {
                        m_internal_k++;
                        m_unsucc_tries = 0;
                }
                return nextPair( m_internal_k );
        }
private:
        std::pair< NodeID, NodeID > nextPair(int k) {
                NodeID lb = k*config.search_space_s;
                NodeID ub = std::min((NodeID) (k+1)*config.search_space_s, m_number_of_nodes);
                std::pair< NodeID, NodeID > ret_value = m_search_space_pointers[k];

                //std::cout <<  "lb " <<  lb <<  " " <<  ub <<  " " <<  m_ub << std::endl;
                if(  m_search_space_pointers[k].second+1 < ub )
                        m_search_space_pointers[k].second++;
                else {
                        if( m_search_space_pointers[k].first + 2 < ub ) {
                                m_search_space_pointers[k].first += 1;
                                m_search_space_pointers[k].second = m_search_space_pointers[k].first + 1;
                        } else {
                                m_search_space_pointers[k].first = lb;
                                m_search_space_pointers[k].second = lb+1;
                        }
                }
                //std::cout <<  ret_value.first <<  " " <<  ret_value.second  << std::endl;

                return ret_value;
        }

        NodeID m_ub;
        NodeID m_internal_k;
        NodeID m_unsucc_tries;
        NodeID m_number_of_nodes;

        PartitionConfig config;
        std::vector< std::pair< NodeID, NodeID > > m_search_space_pointers;
};


#endif /* end of include guard: FULL_SEARCH_SPACE_H49ZQ8A4 */
