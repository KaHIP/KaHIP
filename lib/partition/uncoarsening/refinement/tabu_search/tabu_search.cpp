/******************************************************************************
 * tabu_search.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "data_structure/priority_queues/bucket_pq.h"
#include "quality_metrics.h"
#include "tabu_bucket_queue.h"
#include "tabu_moves_queue.h"
#include "tabu_search.h"
#include "uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement.h"
#include "uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_commons.h"
#include "uncoarsening/refinement/kway_graph_refinement/kway_stop_rule.h"

tabu_search::tabu_search() {
                
}

tabu_search::~tabu_search() {
                
}

EdgeWeight tabu_search::perform_refinement(PartitionConfig & config, graph_access & G, complete_boundary & boundary) {
        quality_metrics qm;
        EdgeWeight input_cut = qm.edge_cut(G);
        EdgeWeight cur_cut   = input_cut;
        EdgeWeight best_cut  = input_cut;
        std::vector< PartitionID > bestmap(G.number_of_nodes(), 0);
        forall_nodes(G, node) {
                bestmap[node] = G.getPartitionIndex(node);
        } endfor
        

        EdgeWeight max_degree        = G.getMaxDegree();
        tabu_bucket_queue* queue     = new tabu_bucket_queue(config, max_degree, G.number_of_nodes());
        tabu_moves_queue* tabu_moves = new tabu_moves_queue();

        matrix* T     = new normal_matrix(G.number_of_nodes(), config.k);
        matrix* gamma = new normal_matrix(G.number_of_nodes(), config.k);

        forall_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        NodeID target            = G.getEdgeTarget(e);
                        PartitionID target_block = G.getPartitionIndex(target);
                        gamma->set_xy( node, target_block, gamma->get_xy(node, target_block) + 1);
                } endfor
        } endfor

        forall_nodes(G, node) {
                bool is_bnd      = false;
                PartitionID pIdx = G.getPartitionIndex(node);
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if( pIdx != G.getPartitionIndex(target)) {
                                is_bnd = true;
                                break;
                        }
                } endfor

                if(is_bnd) {
                        for( unsigned block = 0; block < config.k; block++) {
                                if( gamma->get_xy(node, block) > 0  && G.getPartitionIndex(node) != block) {
                                        queue->insert(node, block, gamma->get_xy(node, block) - gamma->get_xy(node, G.getPartitionIndex(node)));
                                } else {
                                        tabu_moves->insert(node, block, 0);
                                }       
                        }
                }
        } endfor
        
        unsigned no_impro_iterations = 0;
        config.maxT                  = random_functions::nextInt(50, 3000);
        unsigned iteration_limit     = std::min((int)(2*G.number_of_nodes()), 40000);

        std::vector< std::pair<NodeID, PartitionID > > undo_buffer;
        undo_buffer.reserve(G.number_of_edges());

        std::vector<PartitionID> cur_state(G.number_of_nodes());
        int best_idx = -1; int round_counter = -1; unsigned iteration = 0;

        for( iteration = 0, round_counter = 0; iteration < config.maxIter; iteration++) {
                if(!queue->empty()) {
                        Gain gain                          = queue->maxValue();
                        std::pair< NodeID, PartitionID > p = queue->deleteMax();
                        NodeID node                        = p.first;
                        NodeID block                       = p.second;
                        NodeID from                        = G.getPartitionIndex(node);
			
                        if( boundary.getBlockWeight(block) + 1 < config.upper_bound_partition && from != block) {
                                boundary.setBlockWeight(from, boundary.getBlockWeight(from) - 1);
                                boundary.setBlockWeight(block, boundary.getBlockWeight(block) + 1);

                                forall_out_edges(G, e, node) {
                                        NodeID target = G.getEdgeTarget(e);
                                        gamma->set_xy( target, from, gamma->get_xy(target, from) - 1);
                                        gamma->set_xy( target, block, gamma->get_xy(target, block) + 1);
                                } endfor

                                forall_out_edges(G, e, node) {
                                        NodeID target            = G.getEdgeTarget(e);
                                        PartitionID target_block = G.getPartitionIndex(target);

                                        for( unsigned i = 0; i < config.k; i++) {
                                                if(queue->contains( target, i )) {
                                                        if( gamma->get_xy(target, i) == 0) {
                                                                queue->deleteNode(target, i);
                                                        } else {
                                                                queue->changeKey(target, i, gamma->get_xy(target, i) - gamma->get_xy(target, target_block));
                                                        }
                                                } else {
                                                        if( gamma->get_xy(target, i) > 0 && T->get_xy(target, i) < (int)iteration) {
                                                                queue->insert(target, i, gamma->get_xy(target, i) - gamma->get_xy(target, target_block));
                                                        }
                                                }
                                        }
                                } endfor

                                std::pair< NodeID, PartitionID > undo_move;
                                undo_move.first  = node;
                                undo_move.second = G.getPartitionIndex(node);

                                undo_buffer.push_back(undo_move);
                                round_counter++;

                                G.setPartitionIndex(node, block);

                                forall_out_edges(G, e, node) {
                                        for( unsigned i = 0; i < config.k; i++) {
						if(queue->contains( node, i)) {
							if( gamma->get_xy(node, i) ==  0) {
								queue->deleteNode(node,i);
							} else {
								queue->changeKey(node, i, gamma->get_xy(node, i) - gamma->get_xy(node, block));
							}
						} else {
							if(gamma->get_xy(node, i) > 0 && T->get_xy(node, i) < (int)iteration) {
								queue->insert(node, i, gamma->get_xy(node, i) - gamma->get_xy(node, block));
							}
						}
                                        }
                                } endfor

                                cur_cut -= gain;
                        }


                        unsigned tenure = config.maxT;//random_functions::nextInt( config.maxT, 2*config.maxT); 
                        tenure = compute_tenure(iteration, tenure);
                        unsigned small_offset = random_functions::nextInt(1,3);
                        T->set_xy(node, block, iteration + tenure + small_offset);
                        tabu_moves->insert(node, block, iteration + tenure + small_offset);
                        if( T->get_xy( node, from) < (int)iteration ) {
                                T->set_xy(node, from, iteration + tenure);
                                tabu_moves->insert(node, from, iteration + tenure);
                        }

                        if(queue->contains(node, from) ) {
                               queue->deleteNode(node, from);
                        }

                        if(queue->contains(node, block) ) {
                               queue->deleteNode(node, block);
                        }

                }

                //update the best cut found
                if( cur_cut < best_cut ) {
                        best_idx            = undo_buffer.size() - 1 ;
                        best_cut            = cur_cut;
                        no_impro_iterations = 0;
                } else {
                        no_impro_iterations++;
                }


                if( round_counter >= (int)G.number_of_edges() ) {
                      
                       if( best_idx != -1 ) {
                               forall_nodes(G, node) {
                                       cur_state[node] = G.getPartitionIndex(node);
                               } endfor

                               for( int idx = undo_buffer.size()-1; idx > best_idx; idx--) {
                                       G.setPartitionIndex( undo_buffer[idx].first, undo_buffer[idx].second );
                               }
                               forall_nodes(G, node) {
                                       bestmap[node] = G.getPartitionIndex(node); 
                                       G.setPartitionIndex(node, cur_state[node]);
                               } endfor
                       }
                       undo_buffer.clear();
                       best_idx      = -1;
                       round_counter = -1;
                }

                if(no_impro_iterations > iteration_limit) {
                        break;
                }

                //reinsert the buffer
                if( !tabu_moves->empty() ) {
                        while( tabu_moves->minValue() <= (int)iteration ) {
                                std::pair< NodeID, PartitionID > p = tabu_moves->deleteMin();
                                NodeID node = p.first; 
                                NodeID block = p.second; 

                                if( block  == G.getPartitionIndex(node) ) {
                                        unsigned tenure = compute_tenure(iteration, config.maxT);
                                        T->set_xy(node, block, iteration + tenure);
                                        tabu_moves->insert(node, block,iteration + tenure);
                                } else {
					if(gamma->get_xy(node, block) > 0) {
	                                        queue->insert( p.first, p.second,  gamma->get_xy(node, block) - gamma->get_xy(node, G.getPartitionIndex(node)));
					}
                                }
                        }
                }
        }
        if( best_idx != -1 ) {
                forall_nodes(G, node) {
                        cur_state[node] = G.getPartitionIndex(node);
                } endfor


                for( int idx = undo_buffer.size()-1; idx > best_idx; idx--) {
                        G.setPartitionIndex( undo_buffer[idx].first, undo_buffer[idx].second );
                }

                forall_nodes(G, node) {
                        bestmap[node] = G.getPartitionIndex(node); 
                        G.setPartitionIndex(node, cur_state[node]);
                } endfor

        }

        forall_nodes(G, node) {
                G.setPartitionIndex(node, bestmap[node]);
        } endfor

        delete T;
        delete gamma;
        delete queue;
        delete tabu_moves;

        return 0; 
}
