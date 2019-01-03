/******************************************************************************
 * two_way_fm.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "data_structure/priority_queues/bucket_pq.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "macros_assertions.h"
#include "partition_accept_rule.h"
#include "queue_selection_strategie.h"
#include "random_functions.h"
#include "search_stop_rule.h"
#include "tools/quality_metrics.h"
#include "two_way_fm.h"
#include "uncoarsening/refinement/quotient_graph_refinement/partial_boundary.h"

two_way_fm::two_way_fm() {

}

two_way_fm::~two_way_fm() {

}

EdgeWeight two_way_fm::perform_refinement(PartitionConfig & cfg, 
                                          graph_access& G, 
                                          complete_boundary & boundary,
                                          std::vector<NodeID> & lhs_start_nodes, 
                                          std::vector<NodeID> & rhs_start_nodes, 
                                          boundary_pair * pair,        
                                          NodeWeight & lhs_part_weight,
                                          NodeWeight & rhs_part_weight,
                                          EdgeWeight & cut,
                                          bool & something_changed) {

        PartitionConfig config = cfg;//copy it since we make changes on that 
        if(lhs_start_nodes.size() == 0 or rhs_start_nodes.size() == 0) return 0; // nothing to refine

        quality_metrics qm;
        ASSERT_NEQ(pair->lhs, pair->rhs);
        ASSERT_TRUE(assert_directed_boundary_condition(G, boundary, pair->lhs, pair->rhs));
        ASSERT_EQ( cut, qm.edge_cut(G, pair->lhs, pair->rhs));

        refinement_pq* lhs_queue = NULL;
        refinement_pq* rhs_queue = NULL;
        if(config.use_bucket_queues) {
                EdgeWeight max_degree = G.getMaxDegree();
                lhs_queue = new bucket_pq(max_degree); 
                rhs_queue = new bucket_pq(max_degree); 
        } else {
                lhs_queue = new maxNodeHeap(); 
                rhs_queue = new maxNodeHeap(); 
        }

        init_queue_with_boundary(config, G, lhs_start_nodes, lhs_queue, pair->lhs, pair->rhs);  
        init_queue_with_boundary(config, G, rhs_start_nodes, rhs_queue, pair->rhs, pair->lhs);  

        queue_selection_strategy* topgain_queue_select = new queue_selection_topgain(config);
        queue_selection_strategy* diffusion_queue_select = new queue_selection_diffusion(config);
        queue_selection_strategy* diffusion_queue_select_block_target = new queue_selection_diffusion_block_targets(config);
        
        vertex_moved_hashtable moved_idx; 

        std::vector<NodeID> transpositions;

        EdgeWeight inital_cut   = cut;
        int max_number_of_swaps = (int)(boundary.getBlockNoNodes(pair->lhs) + boundary.getBlockNoNodes(pair->rhs));
        int step_limit          = (int)((config.fm_search_limit/100.0)*max_number_of_swaps);
        step_limit              = std::max(step_limit, 15);
        int min_cut_index       = -1;

        refinement_pq* from_queue       = 0;
        refinement_pq* to_queue         = 0;

        PartitionID from                = 0; 
        PartitionID to                  = 0;

        NodeWeight * from_part_weight   = 0;
        NodeWeight * to_part_weight     = 0;

        stop_rule* st_rule                      = new easy_stop_rule();
        partition_accept_rule* accept_partition = NULL;
        if(config.initial_bipartitioning) {
                accept_partition = new ip_partition_accept_rule(config, cut,lhs_part_weight, rhs_part_weight, pair->lhs, pair->rhs);
        } else {
                accept_partition = new normal_partition_accept_rule(config, cut,lhs_part_weight, rhs_part_weight);
        }
        queue_selection_strategy* q_select;

        if(config.softrebalance || config.rebalance || config.initial_bipartitioning) { 
                if(config.initial_bipartitioning) {
                        q_select = diffusion_queue_select_block_target;
                } else {
                        q_select = diffusion_queue_select;
                }
        } else {
                q_select = topgain_queue_select;
        }

        //roll forwards
        EdgeWeight best_cut = cut; 
        int number_of_swaps = 0;
        for(number_of_swaps = 0; number_of_swaps < max_number_of_swaps; number_of_swaps++) {
                if(st_rule->search_should_stop(min_cut_index, number_of_swaps, step_limit)) break;

                if(lhs_queue->empty() && rhs_queue->empty()) { 
                        break; 
                }

                q_select->selectQueue(lhs_part_weight, rhs_part_weight, 
                                pair->lhs, pair->rhs,
                                from,to, 
                                lhs_queue, rhs_queue,
                                &from_queue, &to_queue);

                if(!from_queue->empty()) {
                        Gain gain = from_queue->maxValue();
                        NodeID node = from_queue->deleteMax();

                        ASSERT_TRUE(moved_idx[node].index == NOT_MOVED);

                        boundary.setBlockNoNodes(from, boundary.getBlockNoNodes(from)-1);
                        boundary.setBlockNoNodes(to, boundary.getBlockNoNodes(to)+1);

                        if(from == pair->lhs) {
                                from_part_weight = &lhs_part_weight;         
                                to_part_weight   = &rhs_part_weight;         
                        } else {
                                from_part_weight = &rhs_part_weight;         
                                to_part_weight   = &lhs_part_weight; 
                        }

                        move_node(config, G, node, moved_idx, 
                                        from_queue, to_queue, 
                                        from, to,
                                        pair,
                                        from_part_weight, to_part_weight,
                                        boundary);

                        cut -= gain;

                        if( accept_partition->accept_partition(config, cut, lhs_part_weight, rhs_part_weight, pair->lhs, pair->rhs, config.rebalance)) {
                                ASSERT_TRUE( cut <= best_cut || config.rebalance);
                                if( cut < best_cut ) {
                                        something_changed = true;
                                }
                                best_cut = cut;
                                min_cut_index = number_of_swaps;
                        }

                        transpositions.push_back(node);
                        moved_idx[node].index = MOVED;
                } else {
                        break;
                }

        }

        ASSERT_TRUE(assert_directed_boundary_condition(G, boundary, pair->lhs, pair->rhs)); 
        ASSERT_EQ( cut, qm.edge_cut(G, pair->lhs, pair->rhs));
        
        //roll backwards
        for(number_of_swaps--; number_of_swaps > min_cut_index; number_of_swaps--) {
                ASSERT_TRUE(transpositions.size() > 0);

                NodeID node = transpositions.back();
                transpositions.pop_back();

                PartitionID nodes_partition = G.getPartitionIndex(node);

                if(nodes_partition == pair->lhs) {
                        from_queue       = lhs_queue;
                        to_queue         = rhs_queue;
                        from             = pair->lhs;
                        to               = pair->rhs;
                        from_part_weight = &lhs_part_weight;
                        to_part_weight   = &rhs_part_weight;
                } else {
                        from_queue       = rhs_queue;
                        to_queue         = lhs_queue;
                        from             = pair->rhs;
                        to               = pair->lhs;
                        from_part_weight = &rhs_part_weight;
                        to_part_weight   = &lhs_part_weight;

                }

                boundary.setBlockNoNodes(from, boundary.getBlockNoNodes(from)-1);
                boundary.setBlockNoNodes(to, boundary.getBlockNoNodes(to)+1);

                move_node_back(config, G, node, moved_idx, 
                                from_queue, to_queue, 
                                from, to,
                                pair,
                                from_part_weight, 
                                to_part_weight,
                                boundary);
        }

        //clean up
        cut = best_cut;

        boundary.setEdgeCut(pair, best_cut);
        boundary.setBlockWeight(pair->lhs, lhs_part_weight);
        boundary.setBlockWeight(pair->rhs, rhs_part_weight);

        delete lhs_queue;
        delete rhs_queue;
        delete topgain_queue_select;
        delete diffusion_queue_select;
        delete diffusion_queue_select_block_target;
        delete st_rule;
        delete accept_partition;

        ASSERT_EQ( cut, qm.edge_cut(G, pair->lhs, pair->rhs));
        ASSERT_TRUE(assert_directed_boundary_condition(G, boundary, pair->lhs, pair->rhs)); 
        ASSERT_TRUE(  (int)inital_cut-(int)best_cut >= 0 || cfg.rebalance); 
        // the computed partition shouldnt have a edge cut which is worse than the initial one
        return inital_cut-best_cut;
}


void two_way_fm::move_node(const PartitionConfig & config, 
                           graph_access & G,
                           const NodeID & node,
                           vertex_moved_hashtable & moved_idx,
                           refinement_pq * from_queue,
                           refinement_pq * to_queue,
                           PartitionID from, 
                           PartitionID to,
                           boundary_pair * pair,         
                           NodeWeight * from_part_weight,
                           NodeWeight * to_part_weight,
                           complete_boundary & boundary) {
        //move node
        G.setPartitionIndex(node, to);        
        boundary.deleteNode(node, from, pair);

        EdgeWeight int_degree_node = 0;
        EdgeWeight ext_degree_node = 0;
        bool difficult_update      = int_ext_degree(G, node, to, from, int_degree_node, ext_degree_node);


        if(ext_degree_node > 0) {
                boundary.insert(node, to, pair);
        }

        if(difficult_update)
                boundary.postMovedBoundaryNodeUpdates(node, pair, true, false);


        NodeWeight this_nodes_weight = G.getNodeWeight(node);
        (*from_part_weight) -= this_nodes_weight; 
        (*to_part_weight)   += this_nodes_weight; 

        //update neighbors
        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                PartitionID targets_partition = G.getPartitionIndex(target);

                if((targets_partition != from && targets_partition != to)) {
                        continue; 
                }

                EdgeWeight int_degree = 0;
                EdgeWeight ext_degree = 0;

                PartitionID other_partition = targets_partition == from ? to : from;
                int_ext_degree(G, target, targets_partition, other_partition, int_degree, ext_degree); 

                refinement_pq * queue_to_update = 0;
                if(targets_partition == from) {
                        queue_to_update = from_queue;
                } else {
                        queue_to_update = to_queue;
                }

                Gain gain = ext_degree - int_degree;
                if(queue_to_update->contains(target)) {
                        if(ext_degree == 0) {
                                queue_to_update->deleteNode(target);
                                boundary.deleteNode(target, targets_partition, pair);
                        } else {
                                queue_to_update->changeKey(target, gain);
                        }
                } else {
                        if(ext_degree > 0) {
                                if(moved_idx[target].index == NOT_MOVED) {
                                        queue_to_update->insert(target, gain);
                                }
                                boundary.insert(target, targets_partition, pair);
                        } else {
                                boundary.deleteNode(target, targets_partition, pair);
                        }
                }

        } endfor

}


void two_way_fm::move_node_back(const PartitionConfig & config, 
                                graph_access & G,
                                const NodeID & node,
                                vertex_moved_hashtable & moved_idx,
                                refinement_pq * from_queue,
                                refinement_pq * to_queue,
                                PartitionID from, 
                                PartitionID to,
                                boundary_pair * pair,         
                                NodeWeight * from_part_weight,
                                NodeWeight * to_part_weight,
                                complete_boundary & boundary) {

        ASSERT_NEQ(from, to);
        ASSERT_EQ(from, G.getPartitionIndex(node));

        //move node
        G.setPartitionIndex(node, to);         
        boundary.deleteNode(node, from, pair);

        EdgeWeight int_degree_node = 0;
        EdgeWeight ext_degree_node = 0;
        bool update_difficult = int_ext_degree(G, node, to, from, int_degree_node, ext_degree_node);

        if(ext_degree_node > 0) {
                boundary.insert(node, to, pair);
        }

        if(update_difficult) {
                boundary.postMovedBoundaryNodeUpdates(node, pair, true, false);
        }

        NodeWeight this_nodes_weight = G.getNodeWeight(node);
        (*from_part_weight) -= this_nodes_weight; 
        (*to_part_weight)   += this_nodes_weight; 

        //update neighbors
        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                PartitionID targets_partition = G.getPartitionIndex(target);

                if((targets_partition != from && targets_partition != to)) {
                        //at most difficult update nec.
                        continue; //they dont need to be updated during this refinement
                }

                EdgeWeight int_degree = 0;
                EdgeWeight ext_degree = 0;

                PartitionID other_partition = targets_partition == from ? to : from;
                int_ext_degree(G, target, targets_partition, other_partition, int_degree, ext_degree); 

                if(boundary.contains(target, targets_partition, pair)) {
                        if(ext_degree == 0) {
                                boundary.deleteNode(target, targets_partition, pair);
                        } 
                } else {
                        if(ext_degree > 0) {
                                boundary.insert(target, targets_partition, pair);
                        }
                }

        } endfor

}


void two_way_fm::init_queue_with_boundary(const PartitionConfig & config,
                                          graph_access & G,
                                          std::vector<NodeID> & bnd_nodes,
                                          refinement_pq * queue,                     
                                          PartitionID partition_of_boundary, 
                                          PartitionID other) {

        if(config.permutation_during_refinement == PERMUTATION_QUALITY_FAST) {
                random_functions::permutate_vector_fast(bnd_nodes, false);
        } else if(config.permutation_during_refinement == PERMUTATION_QUALITY_GOOD) {
                random_functions::permutate_vector_good(bnd_nodes, false);
        }


        for( unsigned int i = 0, end = bnd_nodes.size(); i < end; i++) {
                NodeID cur_bnd_node = bnd_nodes[i];
                //compute gain
                EdgeWeight int_degree = 0;
                EdgeWeight ext_degree = 0;

                int_ext_degree(G, cur_bnd_node, partition_of_boundary, other, int_degree, ext_degree);

                Gain gain = ext_degree - int_degree;
                queue->insert(cur_bnd_node, gain);
                ASSERT_TRUE(ext_degree > 0);
                ASSERT_EQ(partition_of_boundary, G.getPartitionIndex(cur_bnd_node));
        }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Assertions for this class////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef NDEBUG 
bool two_way_fm::assert_only_boundary_nodes(graph_access & G, PartialBoundary & lhs_boundary, 
                                            PartitionID lhs, PartitionID rhs) {

        forall_boundary_nodes(lhs_boundary, cur_bnd_node) {
                EdgeWeight int_degree = 0;
                EdgeWeight ext_degree = 0;

                int_ext_degree(G, cur_bnd_node, lhs, rhs, int_degree, ext_degree);

                ASSERT_EQ(G.getPartitionIndex(cur_bnd_node), lhs); 
                ASSERT_TRUE(ext_degree > 0); 
        } endfor
        return true;
}

bool two_way_fm::assert_every_boundary_nodes(graph_access & G, PartialBoundary & lhs_boundary, 
                                             PartitionID lhs, PartitionID rhs) {

        forall_nodes(G, n) {
                EdgeWeight int_degree = 0;
                EdgeWeight ext_degree = 0;
                if(G.getPartitionIndex(n) == lhs) {
                        int_ext_degree(G, n, lhs, rhs, int_degree, ext_degree);

                        if(ext_degree > 0) {
                                ASSERT_TRUE(lhs_boundary.contains(n));
                        }
                }
        } endfor

        return true;
}


bool two_way_fm::assert_directed_boundary_condition(graph_access & G, complete_boundary & boundary, 
                                                    PartitionID lhs, PartitionID rhs) {
        ASSERT_TRUE(assert_only_boundary_nodes(G,  boundary.getDirectedBoundary(lhs, lhs, rhs) , lhs, rhs));
        ASSERT_TRUE(assert_only_boundary_nodes(G,  boundary.getDirectedBoundary(rhs, lhs, rhs) , rhs, lhs));
        ASSERT_TRUE(assert_every_boundary_nodes(G, boundary.getDirectedBoundary(lhs, lhs, rhs) , lhs, rhs));
        ASSERT_TRUE(assert_every_boundary_nodes(G, boundary.getDirectedBoundary(rhs, lhs, rhs) , rhs, lhs));
        return true;
}

#endif
