/******************************************************************************
 * kway_graph_refinement.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <algorithm>
#include <unordered_map>

#include "kway_graph_refinement.h"
#include "kway_graph_refinement_core.h"
#include "kway_stop_rule.h"
#include "quality_metrics.h"
#include "random_functions.h"

kway_graph_refinement::kway_graph_refinement() {
}

kway_graph_refinement::~kway_graph_refinement() {
}

EdgeWeight kway_graph_refinement::perform_refinement(PartitionConfig & config, graph_access & G, 
                                                     complete_boundary & boundary) {

        kway_graph_refinement_core refinement_core;
        
        EdgeWeight overall_improvement = 0;
        int max_number_of_swaps        = (int)(G.number_of_nodes());
        bool sth_changed               = config.no_change_convergence;

        for( unsigned i = 0; i < config.kway_rounds || sth_changed; i++) {
                EdgeWeight improvement = 0;    

                boundary_starting_nodes start_nodes;
                setup_start_nodes(config, G, boundary, start_nodes);

                if(start_nodes.size() == 0) return 0; // nothing to refine

                //metis steplimit
                int step_limit = (int)((config.kway_fm_search_limit/100.0)*max_number_of_swaps);
                step_limit = std::max(step_limit, 15);

                vertex_moved_hashtable moved_idx; 
                improvement += refinement_core.single_kway_refinement_round(config, G, boundary, 
                                                                            start_nodes, step_limit, 
                                                                            moved_idx);

                sth_changed = improvement != 0 && config.no_change_convergence;
                if(improvement == 0) break; 
                overall_improvement += improvement; 

        } 

        ASSERT_TRUE(overall_improvement >= 0); 

        return (EdgeWeight) overall_improvement; 
}

void kway_graph_refinement::setup_start_nodes(PartitionConfig & config, graph_access & G, complete_boundary & boundary,  boundary_starting_nodes & start_nodes) {
        QuotientGraphEdges quotient_graph_edges;
        boundary.getQuotientGraphEdges(quotient_graph_edges);

        std::unordered_map<NodeID, bool> allready_contained;

        for( unsigned i = 0; i < quotient_graph_edges.size(); i++) {
                boundary_pair & ret_value = quotient_graph_edges[i];
                PartitionID lhs           = ret_value.lhs;
                PartitionID rhs           = ret_value.rhs;

                PartialBoundary & partial_boundary_lhs = boundary.getDirectedBoundary(lhs, lhs, rhs);
                forall_boundary_nodes(partial_boundary_lhs, cur_bnd_node) {
                        ASSERT_EQ(G.getPartitionIndex(cur_bnd_node), lhs);
                        if(allready_contained.find(cur_bnd_node) == allready_contained.end() ) { 
                                start_nodes.push_back(cur_bnd_node);
                                allready_contained[cur_bnd_node] = true;
                        }
                } endfor

                PartialBoundary & partial_boundary_rhs = boundary.getDirectedBoundary(rhs, lhs, rhs);
                forall_boundary_nodes(partial_boundary_rhs, cur_bnd_node) {
                        ASSERT_EQ(G.getPartitionIndex(cur_bnd_node), rhs);
                        if(allready_contained.find(cur_bnd_node) == allready_contained.end()) { 
                                start_nodes.push_back(cur_bnd_node);
                                allready_contained[cur_bnd_node] = true;
                        }
                } endfor
        }
}


