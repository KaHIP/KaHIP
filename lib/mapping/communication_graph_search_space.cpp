/******************************************************************************
 * communication_graph_search_space.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "communication_graph_search_space.h"
#include "tools/random_functions.h"

int global_num_nodes = 0;

communication_graph_search_space::communication_graph_search_space(PartitionConfig & config, NodeID number_of_nodes) {
        m_pointer       = 0;
        m_last_pointer  = 0;
        m_unsucc_tries  = 0;
        m_search_deepth = config.communication_neighborhood_dist;
        m_have_to_break = false;
        m_deepth.resize(number_of_nodes);
        global_num_nodes = number_of_nodes;
	this->config = config;
}

void communication_graph_search_space::set_graph_ref( graph_access * C) { 
	this->C = C;
        if( config.communication_neighborhood_dist == 1 ) {
                forall_nodes((*C), node) {
                        forall_out_edges((*C), e, node) {
                                NodeID target = C->getEdgeTarget(e);
                                if( node < target ) {
                                        m_list_of_pairs.push_back( std::pair< NodeID, NodeID> ( node, target ) );
                                }
                        } endfor
                } endfor
        } else {
                std::vector< NodeID > touched_nodes; 
                std::vector< int > deepth(C->number_of_nodes(), -1);
                forall_nodes((*C), node) {
                        touched_nodes.clear();

                        std::queue< NodeID > bfsqueue;
                        bfsqueue.push(node);
                        touched_nodes.push_back(node);

                        deepth[node]  = 0;
                        int cur_deepth  = 0;

                        while(!bfsqueue.empty()) {
                                NodeID v = bfsqueue.front();
                                bfsqueue.pop(); 

                                if (deepth[v] == cur_deepth) {
                                        cur_deepth++;
                                }

                                if( cur_deepth > config.communication_neighborhood_dist ) {
                                        break;
                                }

                                forall_out_edges((*C), e, v) {
                                        NodeID target = C->getEdgeTarget(e);
                                        if(deepth[target] == -1) {
                                                deepth[target] = cur_deepth;
                                                bfsqueue.push(target);
                                                touched_nodes.push_back(target);
                                                m_list_of_pairs.push_back( std::pair< NodeID, NodeID> ( node, target ) );
                                        }
                                } endfor
                        }        

                        for( NodeID v : touched_nodes ) {
                                deepth[v] = -1;
                        }
                } endfor 
        }

        random_functions::permutate_vector_good( m_list_of_pairs);
        m_limit = m_list_of_pairs.size(); 

	for( unsigned int i = 0; i < m_list_of_pairs.size(); i++) {
		m_pair_active[m_list_of_pairs[i]] = true;
	}

        for( unsigned int i = 0; i < m_deepth.size(); i++) {
                m_deepth[i] = -1;
        }


}
communication_graph_search_space::~communication_graph_search_space() {
                
}

