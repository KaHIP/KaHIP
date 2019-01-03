/******************************************************************************
 * quotient_graph_refinement.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <unordered_map>

#include "2way_fm_refinement/two_way_fm.h"
#include "complete_boundary.h"
#include "flow_refinement/two_way_flow_refinement.h"
#include "quality_metrics.h"
#include "quotient_graph_refinement.h"
#include "quotient_graph_scheduling/active_block_quotient_graph_scheduler.h"
#include "quotient_graph_scheduling/simple_quotient_graph_scheduler.h"
#include "uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement.h"
#include "uncoarsening/refinement/kway_graph_refinement/multitry_kway_fm.h"

quotient_graph_refinement::quotient_graph_refinement() {

}

quotient_graph_refinement::~quotient_graph_refinement() {

}

void quotient_graph_refinement::setup_start_nodes(graph_access & G, 
                PartitionID partition, 
                boundary_pair & bp, 
                complete_boundary & boundary,  
                boundary_starting_nodes & start_nodes) {

        start_nodes.resize(boundary.size(partition, &bp));
        NodeID cur_idx = 0;

        PartitionID lhs = bp.lhs;
        PartitionID rhs = bp.rhs;
        PartialBoundary & lhs_b = boundary.getDirectedBoundary(partition, lhs, rhs);

        forall_boundary_nodes(lhs_b, cur_bnd_node) {
                ASSERT_EQ(G.getPartitionIndex(cur_bnd_node), partition);
                start_nodes[cur_idx++] = cur_bnd_node;
        } endfor
}


EdgeWeight quotient_graph_refinement::perform_refinement(PartitionConfig & config, graph_access & G, complete_boundary & boundary) {

        ASSERT_TRUE(boundary.assert_bnodes_in_boundaries());
        ASSERT_TRUE(boundary.assert_boundaries_are_bnodes());

        QuotientGraphEdges qgraph_edges;
        boundary.getQuotientGraphEdges(qgraph_edges);
        quotient_graph_scheduling* scheduler = NULL;

        int factor = ceil(config.bank_account_factor*qgraph_edges.size());
        switch(config.refinement_scheduling_algorithm) {
                case REFINEMENT_SCHEDULING_FAST:
                        scheduler = new simple_quotient_graph_scheduler(config, qgraph_edges, factor);
                        break;
                case REFINEMENT_SCHEDULING_ACTIVE_BLOCKS:
                        scheduler = new active_block_quotient_graph_scheduler(config, qgraph_edges, factor);
                        break;
                case REFINEMENT_SCHEDULING_ACTIVE_BLOCKS_REF_KWAY:
                        scheduler = new active_block_quotient_graph_scheduler(config, qgraph_edges, factor);
                        break;
        }

        EdgeWeight overall_improvement                = 0;
        unsigned int no_of_pairwise_improvement_steps = 0;
        quality_metrics qm;

        do {
                no_of_pairwise_improvement_steps++;
                // ********** preconditions ********************
                ASSERT_TRUE(boundary.assert_bnodes_in_boundaries());
                ASSERT_TRUE(boundary.assert_boundaries_are_bnodes());
                // *************** end *************************

                if(scheduler->hasFinished()) break; //fetch the case where we have no qgraph edges

                boundary_pair & bp = scheduler->getNext();
                PartitionID lhs = bp.lhs;
                PartitionID rhs = bp.rhs;

                NodeWeight lhs_part_weight = boundary.getBlockWeight(lhs);
                NodeWeight rhs_part_weight = boundary.getBlockWeight(rhs);

                EdgeWeight initial_cut_value = boundary.getEdgeCut(&bp);
                if( initial_cut_value < 0 ) continue; // quick fix, for bug 02 (very rare cross combine bug / coarsest level) !

                bool something_changed = false;

#ifndef NDEBUG  
                EdgeWeight oldcut = initial_cut_value;
#endif

                PartitionConfig cfg    = config;
                EdgeWeight improvement = perform_a_two_way_refinement(cfg, G, boundary, bp, 
                                                                      lhs, rhs, 
                                                                      lhs_part_weight, rhs_part_weight, 
                                                                      initial_cut_value, something_changed);

                overall_improvement += improvement;

                EdgeWeight multitry_improvement = 0;
                if(config.refinement_scheduling_algorithm == REFINEMENT_SCHEDULING_ACTIVE_BLOCKS_REF_KWAY ) {
                        multitry_kway_fm kway_ref;
                        std::unordered_map<PartitionID, PartitionID> touched_blocks;

                        multitry_improvement = kway_ref.perform_refinement_around_parts(cfg, G, 
                                                                                boundary, true, 
                                                                                config.local_multitry_fm_alpha, lhs, rhs, 
                                                                                touched_blocks); 

                        if(multitry_improvement > 0) {
                                ((active_block_quotient_graph_scheduler*)scheduler)->activate_blocks(touched_blocks);
                        }

                }

                qgraph_edge_statistics stat(improvement, &bp, something_changed);
                scheduler->pushStatistics(stat);

                //**************** assertions / postconditions ************************** 
                ASSERT_TRUE( oldcut - improvement == qm.edge_cut(G, lhs, rhs) 
                          || config.refinement_scheduling_algorithm == REFINEMENT_SCHEDULING_ACTIVE_BLOCKS_REF_KWAY);
                ASSERT_TRUE(boundary.assert_bnodes_in_boundaries());
                ASSERT_TRUE(boundary.assert_boundaries_are_bnodes());
                ASSERT_TRUE(boundary.getBlockNoNodes(lhs)>0);
                ASSERT_TRUE(boundary.getBlockNoNodes(rhs)>0);
                //*************************** end **************************************** 
        } while(!scheduler->hasFinished());

        delete scheduler;
        return overall_improvement;
}

EdgeWeight quotient_graph_refinement::perform_a_two_way_refinement(PartitionConfig & config, 
                                                                   graph_access & G,
                                                                   complete_boundary & boundary,
                                                                   boundary_pair & bp,
                                                                   PartitionID & lhs, 
                                                                   PartitionID & rhs,
                                                                   NodeWeight & lhs_part_weight,
                                                                   NodeWeight & rhs_part_weight,
                                                                   EdgeWeight & initial_cut_value,
                                                                   bool & something_changed) {

        two_way_fm pair_wise_refinement;
        two_way_flow_refinement pair_wise_flow;

        std::vector<NodeID> lhs_bnd_nodes;
        setup_start_nodes(G, lhs, bp, boundary, lhs_bnd_nodes); 

        std::vector<NodeID> rhs_bnd_nodes; 
        setup_start_nodes(G, rhs, bp, boundary, rhs_bnd_nodes);

        something_changed      = false;
        EdgeWeight improvement = 0;

        quality_metrics qm;
        if(config.refinement_type == REFINEMENT_TYPE_FM_FLOW || config.refinement_type == REFINEMENT_TYPE_FM) {
                improvement = pair_wise_refinement.perform_refinement(config, 
                                                                      G,
                                                                      boundary,
                                                                      lhs_bnd_nodes,
                                                                      rhs_bnd_nodes, 
                                                                      &bp,
                                                                      lhs_part_weight, 
                                                                      rhs_part_weight,
                                                                      initial_cut_value,
                                                                      something_changed); 
                ASSERT_TRUE(improvement >= 0 || config.rebalance); 
        }

        if(config.refinement_type == REFINEMENT_TYPE_FM_FLOW || config.refinement_type == REFINEMENT_TYPE_FLOW){
                lhs_bnd_nodes.clear();
                setup_start_nodes(G, lhs, bp, boundary, lhs_bnd_nodes); 

                rhs_bnd_nodes.clear(); 
                setup_start_nodes(G, rhs, bp, boundary, rhs_bnd_nodes);

                EdgeWeight _improvement = pair_wise_flow.perform_refinement(config, 
                                                                            G, 
                                                                            boundary,
                                                                            lhs_bnd_nodes,
                                                                            rhs_bnd_nodes, 
                                                                            &bp,
                                                                            lhs_part_weight, 
                                                                            rhs_part_weight,
                                                                            initial_cut_value,
                                                                            something_changed);  

                ASSERT_TRUE(_improvement >= 0 || config.rebalance); 
                improvement += _improvement;
        }

        bool only_one_block_is_overloaded = boundary.getBlockWeight(lhs) > config.upper_bound_partition;
        only_one_block_is_overloaded = only_one_block_is_overloaded 
                                     || boundary.getBlockWeight(rhs) > config.upper_bound_partition;
        only_one_block_is_overloaded = only_one_block_is_overloaded && 
                (boundary.getBlockWeight(lhs) <= config.upper_bound_partition || 
                 boundary.getBlockWeight(rhs) <= config.upper_bound_partition);

        if(only_one_block_is_overloaded) {

                PartitionConfig cfg = config;
                cfg.softrebalance   = true;
                cfg.rebalance       = false;

                lhs_bnd_nodes.clear();
                setup_start_nodes(G, lhs, bp, boundary, lhs_bnd_nodes); 

                rhs_bnd_nodes.clear(); 
                setup_start_nodes(G, rhs, bp, boundary, rhs_bnd_nodes);

                improvement += pair_wise_refinement.perform_refinement(cfg, 
                                                                       G,
                                                                       boundary,
                                                                       lhs_bnd_nodes,
                                                                       rhs_bnd_nodes, 
                                                                       &bp,
                                                                       lhs_part_weight, 
                                                                       rhs_part_weight,
                                                                       initial_cut_value,
                                                                       something_changed); 

                ASSERT_TRUE(improvement >= 0 || config.rebalance); 

                if(!config.disable_hard_rebalance && !config.kaffpa_perfectly_balanced_refinement && !config.initial_bipartitioning) {
                                only_one_block_is_overloaded = boundary.getBlockWeight(lhs) > config.upper_bound_partition;
                                only_one_block_is_overloaded = only_one_block_is_overloaded 
                                        || boundary.getBlockWeight(rhs) > config.upper_bound_partition;
                                only_one_block_is_overloaded = only_one_block_is_overloaded && 
                                        (boundary.getBlockWeight(lhs) <= config.upper_bound_partition || 
                                         boundary.getBlockWeight(rhs) <= config.upper_bound_partition);

                        if(only_one_block_is_overloaded) {
                                cfg.softrebalance = true;
                                cfg.rebalance     = true;

                                lhs_bnd_nodes.clear();
                                setup_start_nodes(G, lhs, bp, boundary, lhs_bnd_nodes); 

                                rhs_bnd_nodes.clear(); 
                                setup_start_nodes(G, rhs, bp, boundary, rhs_bnd_nodes);

                                improvement += pair_wise_refinement.perform_refinement(cfg, 
                                                G,
                                                boundary,
                                                lhs_bnd_nodes,
                                                rhs_bnd_nodes, 
                                                &bp,
                                                lhs_part_weight, 
                                                rhs_part_weight,
                                                initial_cut_value,
                                                something_changed); 

                        }
                }
        }

        return improvement;               
}

