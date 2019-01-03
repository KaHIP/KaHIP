/******************************************************************************
 * graph_partitioner.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "coarsening/coarsening.h"
#include "graph_extractor.h"
#include "graph_partitioner.h"
#include "initial_partitioning/initial_partitioning.h"
#include "quality_metrics.h"
#include "tools/random_functions.h"
#include "uncoarsening/uncoarsening.h"
#include "uncoarsening/refinement/mixed_refinement.h"
#include "w_cycles/wcycle_partitioner.h"

graph_partitioner::graph_partitioner() {

}

graph_partitioner::~graph_partitioner() {

}

void graph_partitioner::perform_recursive_partitioning(PartitionConfig & config, graph_access & G) {
        m_global_k = config.k;
        m_global_upper_bound = config.upper_bound_partition;
        m_rnd_bal = random_functions::nextDouble(1,2);
        perform_recursive_partitioning_internal(config, G, 0, config.k-1);
}

void graph_partitioner::perform_recursive_partitioning_internal(PartitionConfig & config, 
                                                                graph_access & G, 
                                                                PartitionID lb, 
                                                                PartitionID ub) {

        G.set_partition_count(2);
        
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // configuration of bipartitioning
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        PartitionConfig bipart_config      = config;
        bipart_config.k                    = 2;
        bipart_config.stop_rule            = STOP_RULE_MULTIPLE_K;
        bipart_config.num_vert_stop_factor = 100;
        double epsilon                     = 0;
        bipart_config.rebalance            = false;
        bipart_config.softrebalance        = true;

        if(config.k < 64) {
                epsilon                     = m_rnd_bal/100.0;
                bipart_config.rebalance     = false;
                bipart_config.softrebalance = false;
        } else {
                epsilon                     = 1/100.0;
        }
        if(m_global_k == 2) {
                epsilon = 3.0/100.0;
        }

        
        bipart_config.upper_bound_partition              = ceil((1+epsilon)*config.work_load/(double)bipart_config.k);
        bipart_config.corner_refinement_enabled          = false;
        bipart_config.quotient_graph_refinement_disabled = false;
        bipart_config.refinement_scheduling_algorithm    = REFINEMENT_SCHEDULING_ACTIVE_BLOCKS;
        bipart_config.kway_adaptive_limits_beta          = log(G.number_of_nodes());
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // end configuration
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        NodeID new_ub_lhs     = floor((lb+ub)/2);
        NodeID new_lb_rhs     = floor((lb+ub)/2+1);
        NodeID num_blocks_lhs = new_ub_lhs - lb + 1;
        NodeID num_blocks_rhs = ub - new_lb_rhs + 1;

        if(config.k % 2 != 0) {
                //otherwise the block weights have to be 
                bipart_config.target_weights.clear();
                bipart_config.target_weights.push_back((1+epsilon)*num_blocks_lhs/(double)(num_blocks_lhs+num_blocks_rhs)*config.work_load);
                bipart_config.target_weights.push_back((1+epsilon)*num_blocks_rhs/(double)(num_blocks_lhs+num_blocks_rhs)*config.work_load);
                bipart_config.initial_bipartitioning  = true;
                bipart_config.refinement_type         = REFINEMENT_TYPE_FM; // flows not supported for odd block weights
        } else {

                bipart_config.target_weights.clear();
                bipart_config.target_weights.push_back(bipart_config.upper_bound_partition);
                bipart_config.target_weights.push_back(bipart_config.upper_bound_partition);
                bipart_config.initial_bipartitioning  = false;
        }

        bipart_config.grow_target = ceil(num_blocks_lhs/(double)(num_blocks_lhs+num_blocks_rhs)*config.work_load);

        perform_partitioning(bipart_config, G);        

        if( config.k > 2 ) {
               graph_extractor extractor;
 
               graph_access extracted_block_lhs;
               graph_access extracted_block_rhs;
               std::vector<NodeID> mapping_extracted_to_G_lhs; // map the new nodes to the nodes in the old graph G
               std::vector<NodeID> mapping_extracted_to_G_rhs; // map the new nodes to the nodes in the old graph G

               NodeWeight weight_lhs_block = 0;
               NodeWeight weight_rhs_block = 0;

               extractor.extract_two_blocks(G, extracted_block_lhs, 
                                               extracted_block_rhs, 
                                               mapping_extracted_to_G_lhs, 
                                               mapping_extracted_to_G_rhs, 
                                               weight_lhs_block, weight_rhs_block);

               PartitionConfig rec_config = config;
               if(num_blocks_lhs > 1) {
                       rec_config.k = num_blocks_lhs;

                       rec_config.largest_graph_weight = weight_lhs_block;
                       rec_config.work_load            = weight_lhs_block;
                       perform_recursive_partitioning_internal( rec_config, extracted_block_lhs, lb, new_ub_lhs);
                       
                       //apply partition
                       forall_nodes(extracted_block_lhs, node) {
                               G.setPartitionIndex(mapping_extracted_to_G_lhs[node], extracted_block_lhs.getPartitionIndex(node));
                       } endfor
 
               } else {
                       //apply partition
                       forall_nodes(extracted_block_lhs, node) {
                               G.setPartitionIndex(mapping_extracted_to_G_lhs[node], lb);
                       } endfor
               }

               if(num_blocks_rhs > 1) {
                       rec_config.k = num_blocks_rhs;
                       rec_config.largest_graph_weight = weight_rhs_block;
                       rec_config.work_load            = weight_rhs_block;
                       perform_recursive_partitioning_internal( rec_config, extracted_block_rhs, new_lb_rhs, ub);

                       forall_nodes(extracted_block_rhs, node) {
                               G.setPartitionIndex(mapping_extracted_to_G_rhs[node], extracted_block_rhs.getPartitionIndex(node));
                       } endfor

               } else {
                       //apply partition
                       forall_nodes(extracted_block_rhs, node) {
                               G.setPartitionIndex(mapping_extracted_to_G_rhs[node], ub);
                       } endfor
               }

        } else {
               forall_nodes(G, node) {
                       if(G.getPartitionIndex(node) == 0) {
                            G.setPartitionIndex(node, lb);
                       } else {
                            G.setPartitionIndex(node, ub);
                       }
               } endfor
        }
       
        G.set_partition_count(config.k);
}

void graph_partitioner::single_run( PartitionConfig & config, graph_access & G) {

        for( unsigned i = 1; i <= config.global_cycle_iterations; i++) {
                PRINT(std::cout <<  "vcycle " << i << " of " << config.global_cycle_iterations  << std::endl;)
                        if(config.use_wcycles || config.use_fullmultigrid)  {
                                wcycle_partitioner w_partitioner;
                                w_partitioner.perform_partitioning(config, G);
                        } else {
                                coarsening coarsen;
                                initial_partitioning init_part;
                                uncoarsening uncoarsen;

                                graph_hierarchy hierarchy;

                                if( config.mode_node_separators ) {
                                        int rnd = random_functions::nextInt(0,3);
                                        if( rnd == 0 ) {
                                                config.edge_rating = SEPARATOR_MULTX;
                                        } else if ( rnd  == 1 ) {
                                                config.edge_rating = WEIGHT;
                                        } else if ( rnd  == 2 ) {
                                                config.edge_rating = SEPARATOR_MAX;
                                        } else if ( rnd  == 3 ) {
                                                config.edge_rating = SEPARATOR_LOG;
                                        } 
                                }
                                coarsen.perform_coarsening(config, G, hierarchy);
                                init_part.perform_initial_partitioning(config, hierarchy);
                                uncoarsen.perform_uncoarsening(config, hierarchy);
                                if( config.mode_node_separators ) {
                                        quality_metrics qm;
                                        std::cout <<  "vcycle result " << qm.separator_weight(G)  << std::endl;
                                }
                        }
                config.graph_allready_partitioned = true;
                config.balance_factor             = 0;
        }
}

void graph_partitioner::perform_partitioning( PartitionConfig & config, graph_access & G) {
        if(config.only_first_level) {
                if( !config.graph_allready_partitioned) {
                        initial_partitioning init_part;
                        init_part.perform_initial_partitioning(config, G);
                }
                
                if( !config.mh_no_mh ) {
                        complete_boundary boundary(&G);
                        boundary.build();
                        refinement* refine      = new mixed_refinement();
                        refine->perform_refinement(config, G, boundary);
                        delete refine;
                }

                return;
        }

        if( config.repetitions == 1 ) {
                single_run(config,G);
        } else {
                quality_metrics qm;
                // currently only for ecosocial
                EdgeWeight best_cut = std::numeric_limits< EdgeWeight >::max();
                std::vector< PartitionID > best_map = std::vector< PartitionID >(G.number_of_nodes());
                for( int i = 0; i < config.repetitions; i++) {
                        forall_nodes(G, node) {
                                G.setPartitionIndex(node,0);
                        } endfor
                        PartitionConfig working_config = config;
                        single_run(working_config, G);

                        EdgeWeight cur_cut = qm.edge_cut(G);
                        if( cur_cut < best_cut ) {
                                forall_nodes(G, node) {
                                        best_map[node] = G.getPartitionIndex(node);
                                } endfor

                                best_cut = cur_cut;
                        }
                }

                forall_nodes(G, node) {
                        G.setPartitionIndex(node, best_map[node]);
                } endfor

        }
}

