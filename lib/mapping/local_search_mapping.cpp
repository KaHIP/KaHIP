/******************************************************************************
 * local_search_mapping.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "tools/quality_metrics.h"
#include "tools/random_functions.h"
#include "tools/timer.h"
#include "local_search_mapping.h"
#include "full_search_space.h"
#include "communication_graph_search_space.h"

local_search_mapping::local_search_mapping() {

}

local_search_mapping::~local_search_mapping() {

}


bool local_search_mapping::perform_single_swap(graph_access & C, matrix & D, std::vector< NodeID > & perm_rank, NodeID swap_lhs, NodeID swap_rhs) {
        NodeWeight old_volume      = total_volume;
        NodeWeight old_lhs_contrib = node_contribution[swap_lhs];
        NodeWeight old_rhs_contrib = node_contribution[swap_rhs];

        // we multiply by two since contributions are on both sides
        total_volume -= 2*node_contribution[swap_lhs];
        total_volume -= 2*node_contribution[swap_rhs];

        // fix adjacent candiates
        forall_out_edges(C, e, swap_lhs) {
                NodeID target = C.getEdgeTarget(e);
                if( target == swap_rhs ) {
                        NodeWeight comm_vol     = C.getEdgeWeight(e);
                        NodeID perm_rank_node   = perm_rank[swap_lhs];
                        NodeID perm_rank_target = perm_rank[swap_rhs];
                        NodeWeight cur_vol      = comm_vol*D.get_xy(perm_rank_node, perm_rank_target);
                        total_volume += 2*cur_vol;
                        break;
                }
        } endfor

        node_contribution[swap_lhs] = 0;
        node_contribution[swap_rhs] = 0;

        std::swap(perm_rank[swap_lhs], perm_rank[swap_rhs]);
        update_node_contribution( C, D, perm_rank, swap_lhs, swap_rhs );

        total_volume += 2*node_contribution[swap_lhs]; 
        total_volume += 2*node_contribution[swap_rhs]; 

        // fix adjacent candiates
        forall_out_edges(C, e, swap_lhs) {
                NodeID target = C.getEdgeTarget(e);
                if( target == swap_rhs ) {
                        NodeWeight comm_vol     = C.getEdgeWeight(e);
                        NodeID perm_rank_node   = perm_rank[swap_lhs];
                        NodeID perm_rank_target = perm_rank[swap_rhs];
                        NodeWeight cur_vol      = comm_vol*D.get_xy(perm_rank_node, perm_rank_target);
                        total_volume -= 2*cur_vol;
                        break;
                }
        } endfor

        if( total_volume < old_volume ) {
                PRINT(std::cout <<  "log> improvement " <<  total_volume <<  " " <<  old_volume << std::endl;)
                return true;
        } else {
                std::swap(perm_rank[swap_lhs], perm_rank[swap_rhs]);
                update_node_contribution( C, D, perm_rank, swap_lhs, swap_rhs );
                node_contribution[swap_lhs] = old_lhs_contrib;
                node_contribution[swap_rhs] = old_rhs_contrib;
                total_volume = old_volume;
                return false;
        }
}

void local_search_mapping::update_node_contribution( graph_access & C, matrix & D, std::vector< NodeID > & perm_rank, NodeID swap_lhs, NodeID swap_rhs) {
        forall_out_edges(C, e, swap_lhs) {
                NodeID target                   = C.getEdgeTarget(e);
                NodeWeight comm_vol             = C.getEdgeWeight(e);
                NodeID perm_rank_node           = perm_rank[swap_lhs];
                NodeID perm_rank_target         = perm_rank[target];
                NodeWeight cur_vol              = comm_vol*D.get_xy(perm_rank_node, perm_rank_target);
                node_contribution[ swap_lhs ]  += cur_vol;

                // update adjacent node contrib
                if( target != swap_rhs) {
                        node_contribution[ target ] -= comm_vol*D.get_xy(perm_rank[swap_rhs], perm_rank_target);
                        node_contribution[ target ] += cur_vol;
                }
        } endfor
        forall_out_edges(C, e, swap_rhs) {
                NodeID target                   = C.getEdgeTarget(e);
                NodeWeight comm_vol             = C.getEdgeWeight(e);
                NodeID perm_rank_node           = perm_rank[swap_rhs];
                NodeID perm_rank_target         = perm_rank[target];
                NodeWeight cur_vol              = comm_vol*D.get_xy(perm_rank_node, perm_rank_target);
                node_contribution[ swap_rhs ]  += cur_vol;

                if( target != swap_lhs) {
                        node_contribution[ target ] -= comm_vol*D.get_xy(perm_rank[swap_lhs], perm_rank_target);
                        node_contribution[ target ] += cur_vol;
                }
        } endfor
}
