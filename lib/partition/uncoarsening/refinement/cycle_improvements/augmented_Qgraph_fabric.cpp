/******************************************************************************
 * augmented_Qgraph_fabric.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <algorithm>

#include "algorithms/cycle_search.h"
#include "augmented_Qgraph_fabric.h"
#include "data_structure/priority_queues/bucket_pq.h"
#include "partition_snapshooter.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_core.h"
#include "uncoarsening/refinement/kway_graph_refinement/kway_stop_rule.h"
#include "uncoarsening/refinement/quotient_graph_refinement/2way_fm_refinement/vertex_moved_hashtable.h"

augmented_Qgraph_fabric::augmented_Qgraph_fabric() {
}

augmented_Qgraph_fabric::~augmented_Qgraph_fabric() {

}

void augmented_Qgraph_fabric::cleanup_eligible() {
        for( unsigned i = 0; i < m_tomake_eligible.size(); i++) {
                m_eligible[m_tomake_eligible[i]] = true;
        }
        m_tomake_eligible.clear();
 }

bool augmented_Qgraph_fabric::build_augmented_quotient_graph( PartitionConfig & config, 
                                                              graph_access & G, 
                                                              complete_boundary & boundary, 
                                                              augmented_Qgraph & aqg, 
                                                              unsigned & s, bool rebalance, bool plus) {

        graph_access G_bar;
        boundary.getUnderlyingQuotientGraph(G_bar); 
        if(m_eligible.size() != G.number_of_nodes()) {
                m_eligible.resize(G.number_of_nodes());
                forall_nodes(G, node) {
                        m_eligible[node] = true;
                } endfor
        } else {
                cleanup_eligible();
        }

        if(!rebalance) {
                std::vector<block_pair_difference> vec_bpd;
                forall_nodes(G_bar, lhs) {
                        forall_out_edges(G_bar, e, lhs) {
                                EdgeID rhs = G_bar.getEdgeTarget(e);

                                block_pair_difference bpd;
                                bpd.lhs = lhs;
                                bpd.rhs = rhs;
                                vec_bpd.push_back(bpd);
                        } endfor 
                } endfor


                for( unsigned j = 0; j < config.kaba_packing_iterations; j++) {
                        random_functions::permutate_vector_good_small(vec_bpd);
                        bool variant_to_use = plus;
                        for( unsigned i = 0; i < vec_bpd.size(); i++) {
                                boundary_pair bp;
                                bp.k   = config.k;
                                bp.lhs = vec_bpd[i].lhs;
                                bp.rhs = vec_bpd[i].rhs;

                                if( plus && config.kaba_flip_packings) {
                                        //best of both worlds
                                        variant_to_use = random_functions::nextBool();
                                }

                                local_search( config, variant_to_use, G, boundary, aqg, bp, s);
                        }
                }
        } else {
                std::vector<block_pair_difference> vec_bpd;
                bool graph_model_will_be_feasable = false;
                forall_nodes(G_bar, lhs) {
                        forall_out_edges(G_bar, e, lhs) {
                                EdgeID rhs = G_bar.getEdgeTarget(e);

                                block_pair_difference bpd;
                                bpd.lhs = lhs;
                                bpd.rhs = rhs;
                                vec_bpd.push_back(bpd);

                                //make the underlying model feasable
                                if( boundary.getBlockWeight(lhs) > config.upper_bound_partition 
                                 && boundary.getBlockWeight(rhs) < config.upper_bound_partition 
                                 && !graph_model_will_be_feasable) {
                                        boundary_pair bp;
                                        bp.k   = config.k;
                                        bp.lhs = lhs;
                                        bp.rhs = rhs;

                                        bool success = local_search( config, false, G, boundary, aqg, bp, s);

                                        if( success ) {
                                                graph_model_will_be_feasable = true;
                                        } // the else case can happen if the quotient graph data structure is not up to date 
                                }
                        } endfor 
                } endfor

                if( !graph_model_will_be_feasable) {
                        // fall back solution
			std::deque<NodeID>* bfsqueue = new std::deque<NodeID>;
                        std::vector< int > parent(G_bar.number_of_nodes(), -1); 

                        std::vector<NodeID> start_vertices;
                        std::vector<NodeID> candidates;
                        forall_nodes(G_bar, lhs) {
                                if( boundary.getBlockWeight(lhs) > config.upper_bound_partition ) {
                                        start_vertices.push_back(lhs);
                                } else if ( boundary.getBlockWeight(lhs) < config.upper_bound_partition) {
                                        candidates.push_back(lhs);
                                }
                        } endfor

                        random_functions::permutate_vector_good_small(start_vertices);
                        for( unsigned i = 0; i < start_vertices.size(); i++) {
                                bfsqueue->push_back(start_vertices[i]);
                                parent[start_vertices[i]] = start_vertices[i];
                        }

                        while(!bfsqueue->empty()) {
                                NodeID lhs = bfsqueue->front();
                                bfsqueue->pop_front();

                                forall_out_edges(G_bar, e, lhs) {
                                        NodeID rhs = G_bar.getEdgeTarget(e);

                                        if(parent[rhs] == -1 && boundary.getDirectedBoundary(lhs, lhs, rhs).size() > 0) {
                                                parent[rhs] = lhs;
                                                bfsqueue->push_back(rhs);
                                        }
                                } endfor
                        }

                        delete bfsqueue;
                    
                        int cur_block;
                        int start_block;
                        bool candiate_set_was_empty = false;
                        std::vector<NodeID> tmp_candidates;
                        tmp_candidates = candidates; // for the connected component case
                        do {
                                if(candidates.size() == 0) {
                                        candiate_set_was_empty = true;
                                        break;
                                }
                                unsigned int r_idx = random_functions::nextInt(0, candidates.size()-1);
                                cur_block          = candidates[r_idx];
                                std::swap(candidates[r_idx], candidates[candidates.size()-1]);
                                candidates.pop_back(); 

                        } while( parent[cur_block] == -1 ); // in this case the vertex is not reachable 
                        
                        //special case for more connected components, 
                        //move a random node from an overloaded block to cur_block (which is a connected component)
                        if(candiate_set_was_empty) {
                                unsigned int r_idx    = random_functions::nextInt(0, tmp_candidates.size()-1);
                                PartitionID cur_block = tmp_candidates[r_idx];

                                do {
                                        unsigned int node       = random_functions::nextInt(0, G.number_of_nodes()-1);
                                        PartitionID nodes_block = G.getPartitionIndex(node);
                                        if( nodes_block != cur_block 
                                         && boundary.getBlockWeight(nodes_block) > config.upper_bound_partition) {
                                               PartitionID from = G.getPartitionIndex(node);
                                               PartitionID to   = cur_block;
                                               perform_simple_move( config, G, boundary, node,from, to);
                                               return true;
                                        }
                                } while( true );
                        } // else

                        start_block = cur_block;
                        std::unordered_map< PartitionID, bool > allready_performed_local_search;

                        while( boundary.getBlockWeight( cur_block ) <= config.upper_bound_partition ) {
                                boundary_pair bp;
                                bp.k   = config.k;
                                bp.lhs = parent[cur_block];
                                bp.rhs = cur_block;

                                cur_block    = parent[cur_block];
                                bool success = local_search( config, false, G, boundary, aqg, bp, 1);

                                if(!success) {
                                        candidates.push_back( start_block );
                                        rebalance_fall_back(config, G, G_bar, boundary, candidates, parent, aqg); // this happens on dense graphs if no node was eligible
                                        return true;
                                } 
                                allready_performed_local_search[config.k*bp.lhs+bp.rhs] = true;
                        } 
                        for( unsigned j = 0; j < config.kaba_packing_iterations; j++) {
                                random_functions::permutate_vector_good_small(vec_bpd);

                                for( unsigned i = 0; i < vec_bpd.size(); i++) {
                                        boundary_pair bp;
                                        bp.k   = config.k;
                                        bp.lhs = vec_bpd[i].lhs;
                                        bp.rhs = vec_bpd[i].rhs;

                                        local_search( config, false, G, boundary, aqg, bp, s);
                                }
                        }
                } else {
                        for( unsigned j = 0; j < config.kaba_packing_iterations; j++) {
                                random_functions::permutate_vector_good_small(vec_bpd);

                                for( unsigned i = 0; i < vec_bpd.size(); i++) {
                                        boundary_pair bp;
                                        bp.k   = config.k;
                                        bp.lhs = vec_bpd[i].lhs;
                                        bp.rhs = vec_bpd[i].rhs;

                                        local_search( config, false, G, boundary, aqg, bp, s);
                                }
                        }
                }

        }

        return false;
}

bool augmented_Qgraph_fabric::construct_local_searches_on_qgraph_edge( PartitionConfig & config, graph_access & G, 
                                                                       complete_boundary & boundary, augmented_Qgraph & aqg, 
                                                                       boundary_pair & pair, 
                                                                       unsigned s,
                                                                       bool plus) {
        PartitionID lhs = pair.lhs;
        PartitionID rhs = pair.rhs;

        //initialize todo list
        PartialBoundary & lhs_b = boundary.getDirectedBoundary(lhs, lhs, rhs);
        std::vector<NodeID> lhs_boundary; // todo list
        forall_boundary_nodes(lhs_b, node) {
                if(m_eligible[node]) {
                        lhs_boundary.push_back(node);
                }
        } endfor

        if(lhs_boundary.size() == 0) {  
                //nothing todo 
                return false; 
        }

        //commons = kway_graph_refinement_commons::getInstance(config);
        for( unsigned i = 0; i < 1; i++) {

                if(lhs_boundary.size() == 0) return false;

                pairwise_local_search pls;

                NodeID start_node = lhs_boundary[0];
                find_eligible_start_node( G, lhs, rhs,  lhs_boundary, m_eligible, start_node);
                
                if(!m_eligible[start_node]) return false; // in this case the lhs_boundary was empty and we cant move a node

                if(plus) {
                        more_locallized_search(config, G,  boundary, lhs, rhs, start_node, s, pls);
                } else {
                        directed_more_locallized_search(config, G,  boundary, lhs, rhs, start_node, s, pls);
                }

                aqg.commit_pairwise_local_search(pair, pls);

                if( plus ) {
                        // keep things simple
                        boundary_pair opp_pair = pair;
                        std::swap(opp_pair.lhs, opp_pair.rhs);
                        aqg.commit_pairwise_local_search(opp_pair, pls);
                }
        }
        return true;
}

//this method performes a directed localized local search and UNDOs them
//these searches are for the augmented qgraph structure for balanced graph partitioning
void augmented_Qgraph_fabric::directed_more_locallized_search(PartitionConfig & config, graph_access & G, 
                                                   complete_boundary & boundary,  
                                                   PartitionID & lhs, PartitionID & rhs,
                                                   NodeID start_node, unsigned & number_of_swaps, pairwise_local_search & pls) {

        //commons = kway_graph_refinement_commons::getInstance(config);
        EdgeWeight max_degree = G.getMaxDegree();
        refinement_pq* queue  = new bucket_pq(max_degree);

        EdgeWeight int_degree = 0;
        EdgeWeight ext_degree = 0;
        m_twfm.int_ext_degree(G, start_node, lhs, rhs, int_degree, ext_degree);

        Gain gain = ext_degree - int_degree; 
        queue->insert(start_node, gain);
        
        if(queue->empty()) {delete queue; return;}

        ////roll forwards
        int movements    = 0;
        int overall_gain = 0;

        kway_stop_rule* stopping_rule = new kway_simple_stop_rule(config);

        int min_cut_index    = 0;
        int step_limit       = 200;
        EdgeWeight input_cut = boundary.getEdgeCut(lhs, rhs);
        EdgeWeight min_cut   = input_cut;
        PartitionID from     = lhs;
        PartitionID to       = rhs;

        for(movements = 0; movements < (int)number_of_swaps; movements++) {
                if( queue->empty() ) {
                        break;
                }
                if( stopping_rule->search_should_stop(min_cut_index, movements, step_limit) ) break;


                Gain gain   = queue->maxValue();
                NodeID node = queue->deleteMax();

                move_node(config, G, node, queue, boundary, from, to);

                overall_gain += gain;
                input_cut -= gain;
        
                stopping_rule->push_statistics(gain);

                if(input_cut < min_cut) {
                        min_cut_index = movements;
                        min_cut       = input_cut;
                }

                pls.vertex_movements.push_back(node);
                pls.block_movements.push_back(to);
                pls.gains.push_back(overall_gain);
                m_tomake_eligible.push_back(node);
        }

        ////roll backwards
        int idx = pls.vertex_movements.size()-1;
        for(; idx >= 0; idx--) {
                NodeID node = pls.vertex_movements[idx];
                move_node(config, G, node,  queue, boundary, to, from);

                //block the neighboring nodes to avoid conflicts
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if(m_eligible[target]) m_tomake_eligible.push_back(target);
                        m_eligible[target] = false;
                } endfor
        }
        delete queue;
        delete stopping_rule;
}

//this method performes a directed localized local search and UNDOs them
//these searches are for the augmented qgraph structure for balanced graph partitioning
void augmented_Qgraph_fabric::more_locallized_search(PartitionConfig & config, graph_access & G, 
                                                   complete_boundary & boundary,
                                                   PartitionID & lhs, PartitionID & rhs,
                                                   NodeID start_node, unsigned & number_of_swaps, pairwise_local_search & pls) {

        //commons = kway_graph_refinement_commons::getInstance(config);
        refinement_pq* queue_lhs = NULL;
        refinement_pq* queue_rhs = NULL;
        EdgeWeight max_degree = G.getMaxDegree();
        queue_lhs = new bucket_pq(max_degree);
        queue_rhs = new bucket_pq(max_degree);

        EdgeWeight int_degree = 0;
        EdgeWeight ext_degree = 0;
        m_twfm.int_ext_degree(G, start_node, lhs, rhs, int_degree, ext_degree);

        Gain gain = ext_degree - int_degree; 
        queue_lhs->insert(start_node, gain);

        //=====================================
        // find a start node for the rhs queue
        //=====================================
        NodeID start_node_rhs = start_node; // some dummy initialization
        Gain   max_gain = std::numeric_limits<Gain>::min();
        forall_out_edges(G, e, start_node) {
                NodeID target = G.getEdgeTarget(e);
                if( G.getPartitionIndex(target) == rhs && m_eligible[target]) {
                        m_twfm.int_ext_degree(G, target, rhs, lhs, int_degree, ext_degree);
                        if( ext_degree - int_degree > max_gain ) {
                                max_gain = ext_degree - int_degree;
                                start_node_rhs = target;
                        } 
                }
        } endfor
        //=====================================

        if( m_eligible[start_node_rhs] && start_node_rhs != start_node) {
                queue_rhs->insert(start_node_rhs, max_gain);
        }
        
        if(queue_lhs->empty() || queue_rhs->empty()) {delete queue_lhs; delete queue_rhs; return;}
        // queues initalized

        ////roll forwards
        int movements = 0;
        int overall_gain = 0;

        kway_stop_rule* stopping_rule = new kway_simple_stop_rule(config);

        int min_cut_index    = 0;
        int step_limit       = 200;
        EdgeWeight input_cut = boundary.getEdgeCut(lhs, rhs);
        EdgeWeight min_cut   = input_cut;
        PartitionID from     = lhs;
        PartitionID to       = rhs;

        refinement_pq * queue = NULL;
        refinement_pq * to_queue = NULL;

        int diff = 0;

        for(movements = 0; movements < (int)number_of_swaps; movements++) {
                if( queue_lhs->empty() || queue_rhs->empty()) {
                        break;
                }
                if( stopping_rule->search_should_stop(min_cut_index, movements, step_limit) ) break;


                Gain gain_lhs   = queue_lhs->maxValue();
                Gain gain_rhs   = queue_rhs->maxValue();

                bool coin = false;
                switch(config.kaba_lsearch_p) {
                case COIN_DIFFTIE:
                        coin = random_functions::nextBool();
                        if(coin) {
                                if( gain_rhs > gain_lhs ) {
                                        queue = queue_rhs;
                                } else if ( gain_rhs < gain_lhs ) {
                                        queue = queue_lhs;
                                } else {
                                        queue = queue_lhs;
                                }
                        } else {
                                queue = queue_lhs;
                        }
                        break;

                case COIN_RNDTIE:
                        coin = random_functions::nextBool();
                        if(coin) {
                                if( gain_rhs > gain_lhs ) {
                                        queue = queue_rhs;
                                } else if ( gain_rhs < gain_lhs ) {
                                        queue = queue_lhs;
                                } else {
                                        coin = random_functions::nextBool();
                                        if(coin) {
                                                queue = queue_rhs;
                                        } else {
                                                queue = queue_lhs;
                                        }
                                }
                        } else {
                                queue = queue_lhs;
                        }
                        break;
                case NOCOIN_DIFFTIE:
                        if( gain_rhs > gain_lhs ) {
                                queue = queue_rhs;
                        } else if ( gain_rhs < gain_lhs ) {
                                queue = queue_lhs;
                        } else {
                                queue = queue_lhs;
                        }
                        break;

                case NOCOIN_RNDTIE:
                        if( gain_rhs > gain_lhs ) {
                                queue = queue_rhs;
                        } else if ( gain_rhs < gain_lhs ) {
                                queue = queue_lhs;
                        } else {
                                coin = random_functions::nextBool();
                                if(coin) {
                                        queue = queue_rhs;
                                } else {
                                        queue = queue_lhs;
                                }
                        }
                        break;
                }

                NodeID node = queue->deleteMax();
                if( queue == queue_rhs ) {
                        from     = rhs;
                        to       = lhs;
                        gain     = gain_rhs;
                        to_queue = queue_lhs;
                        diff    += G.getNodeWeight(node);
                } else {
                        from     = lhs;
                        to       = rhs;
                        gain     = gain_lhs;
                        to_queue = queue_rhs;
                        diff    -= G.getNodeWeight(node);
                }


                move_node(config, G, node,  queue, to_queue, boundary, from, to);

                overall_gain += gain;
                input_cut    -= gain;
        
                stopping_rule->push_statistics(gain);

                if(input_cut < min_cut && diff == 0) {
                        min_cut = input_cut;

                        pls.vertex_movements.clear();
                        pls.block_movements.clear();
                        pls.gains.clear();
                } else {
                        pls.vertex_movements.push_back(node);
                        pls.block_movements.push_back(to);
                        pls.gains.push_back(overall_gain);
                }
                m_tomake_eligible.push_back(node);
        }

        ////roll backwards
        int idx = pls.vertex_movements.size()-1;
        //int idx = movements-1;
        for(; idx >= 0; idx--) {
                NodeID node      = pls.vertex_movements[idx];
                PartitionID from = G.getPartitionIndex(node);
                PartitionID to   = from == lhs ? rhs : lhs;
                perform_simple_move( config, G, boundary, node, from, to);

                //block the neighboring nodes to avoid conflicts
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if(m_eligible[target]) m_tomake_eligible.push_back(target);
                        m_eligible[target] = false;
                } endfor
        }

        delete queue_lhs;
        delete queue_rhs;
        delete stopping_rule;
}

