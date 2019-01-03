/******************************************************************************
 * advanced_models.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef ADVANCED_MODELS_PR6SXN3G
#define ADVANCED_MODELS_PR6SXN3G

#include <vector>

#include "augmented_Qgraph_fabric.h"
#include "definitions.h"
#include "quality_metrics.h"
#include "uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_commons.h"
#include "uncoarsening/refinement/quotient_graph_refinement/2way_fm_refinement/two_way_fm.h"
#include "uncoarsening/refinement/refinement.h"

class advanced_models {
        public:
                advanced_models();
                virtual ~advanced_models();

                bool compute_vertex_movements_rebalance( PartitionConfig & config, 
                                graph_access & G, 
                                complete_boundary & boundary, 
                                augmented_Qgraph & aqg, 
                                unsigned & s);

                bool compute_vertex_movements_rebalance_ultra( PartitionConfig & config, 
                                graph_access & G, 
                                complete_boundary & boundary, 
                                augmented_Qgraph & aqg, 
                                unsigned & s);

                bool compute_vertex_movements_ultra_model( PartitionConfig & config, 
                                graph_access & G, 
                                complete_boundary & boundary, 
                                augmented_Qgraph & aqg, 
                                unsigned & s, bool zero_weight_cycle);

                void perform_augmented_move( PartitionConfig & config, 
                                graph_access & G, 
                                complete_boundary & boundary, 
                                std::vector<NodeID> & cycle, 
                                NodeID & s, NodeID & t, 
                                augmented_Qgraph & aqg);

                static unsigned long conflicts;
        private:
                inline
                        bool build_ultra_model( PartitionConfig & config, 
                                        graph_access & G, 
                                        graph_access & G_bar,
                                        complete_boundary & boundary, 
                                        augmented_Qgraph & aqg, 
                                        std::vector<bool> & feasable_edge,
                                        unsigned & steps,
                                        graph_access & cycle_problem, NodeID & s,
                                        std::vector<NodeID> & id_mapping);

                inline
                        bool build_rebalance_model( PartitionConfig & config, 
                                        graph_access & G, 
                                        graph_access & G_bar,
                                        complete_boundary & boundary, 
                                        augmented_Qgraph & aqg, 
                                        std::vector<bool> & feasable_edge,
                                        unsigned & steps,
                                        graph_access & cycle_problem, NodeID & s, NodeID & t,
                                        std::vector<NodeID> & id_mapping);


                inline 
                        bool handle_ultra_model_conflicts( PartitionConfig & config, 
                                        graph_access & cycle_problem, 
                                        complete_boundary & boundary,
                                        std::vector<NodeID> & id_mapping, 
                                        std::vector<bool> & feasable_edge, 
                                        std::vector< NodeID > & cycle, 
                                        NodeID & s, augmented_Qgraph & aqg, bool remove_only_between_layers = false);

                inline 
                        bool cycleorpath_has_conflicts( PartitionConfig & config, 
                                        complete_boundary & boundary,
                                        std::vector< NodeID > & cycleorpath, 
                                        NodeID & s, augmented_Qgraph & aqg);



};

inline 
bool advanced_models::cycleorpath_has_conflicts( PartitionConfig & config, 
                complete_boundary & boundary,
                std::vector< NodeID > & cycleorpath, 
                NodeID & s, augmented_Qgraph & aqg) {
        //check_block_overload
        //check_dependant_edges (depandant local searches)
        //check_augmented_circle_is_simple
        //check_block_get_empty
        // ============ detect wether we have a conflict
        std::vector<bool> block_seen( config.k, false );
        bool conflict_detected = false;
        for( unsigned i = 0; i < cycleorpath.size()-1; i++) {
                if(cycleorpath[i] == s ) {
                        continue;
                }
                if( cycleorpath[i] == cycleorpath[i+1] + config.k || cycleorpath[i]+config.k == cycleorpath[i+1]) {
                        //augmented model case 
                        continue;
                }

                if( cycleorpath[i+1] == s ) {
                        if( cycleorpath.size() > 3 ) { // we dont have a cycle s -> block -> s 
                                // the edge we are interested in is the edge R; it decides how much vertices cycle[i] will take
                                // cycle[i-1] ->_R  cycle[i] -> s
                                bool pos_found = false;
                                unsigned rhs_pos = i;
                                unsigned cur_pos = i;
                                div_t rhs_tmp = div((int)cycleorpath[i],(int)config.k);
                                unsigned num_vert_for_rhs;

                                do {
                                        unsigned lhs_pos = cur_pos >= 1 ? cur_pos - 1 : cycleorpath.size() - 2; 
                                        if(lhs_pos == rhs_pos) break; // not a "real" circle (no edge in the circle is associated with vertex movements)

                                        div_t lhs_tmp = div((int)cycleorpath[lhs_pos],(int)config.k);
                                        if( lhs_tmp.rem != rhs_tmp.rem) {
                                                pos_found = true;
                                                num_vert_for_rhs = lhs_tmp.quot+1;
                                                break;
                                        } else {
                                                if( cur_pos != 0 ) {
                                                        cur_pos--;
                                                } else {
                                                        cur_pos = cycleorpath.size();
                                                }
                                        }
                                        
                                } while(true);
                                if(!pos_found) continue;

                                if( !(boundary.getBlockWeight( rhs_tmp.rem ) + num_vert_for_rhs <= config.upper_bound_partition)) {
                                        conflict_detected = true;
                                        break;
                                } 
                        }
                }

                div_t lhs_tmp = div((int)cycleorpath[i],(int)config.k);
                PartitionID block = lhs_tmp.rem;

                if( block_seen[ block ] ) {
                        conflict_detected = true;
                        break;
                } else {
                        block_seen[block] = true;
                }
        } //postcondition: each block touched by the cycle has block_seen == true

        if(!conflict_detected) {
                //check for cycles of length two that might contain a conflict
                //check wether s is contained
                bool s_contained = false;
                for( unsigned i = 0; i < cycleorpath.size(); i++) {
                        if( cycleorpath[i] == s ) {
                                s_contained = true;
                                break;

                        }
                }

                if( ! s_contained ) { 
                        PartitionID lhs    = 0, rhs       = 0;
                        unsigned counter   = 0;

                        for( unsigned i = 0; i < config.k; i++) {
                                if(block_seen[i]) {
                                        if( counter == 0 ) {
                                                lhs = i;
                                        } else if ( counter == 1 ) {
                                                rhs = i;
                                        }
                                        counter++;
                                }
                                if(counter > 2) break;
                        }

                        if(counter == 2) {
                                // find the right layers
                                unsigned layer_lhs = 0, layer_rhs = 0;
                                for( unsigned i = 0; i < cycleorpath.size(); i++) {
                                        div_t tmp = div((int)cycleorpath[i],(int)config.k);
                                        if( tmp.rem == (int)lhs ) {
                                                if( tmp.quot > (int)layer_lhs) {
                                                        layer_lhs = tmp.quot;
                                                }
                                        } else if (tmp.rem == (int)rhs ) {
                                                if( tmp.quot > (int)layer_rhs) {
                                                        layer_rhs = tmp.quot;
                                                }
                                        }
                                }
                                conflict_detected = aqg.check_conflict( config, lhs, rhs, layer_lhs + 1, layer_rhs + 1);
                        }
                } // otherwise we have a path; a path cannot be conflicted in this case

        }
        if(!conflict_detected) {
                //check wether a block would get empty
                unsigned idx = 0;
                for( unsigned i = 0; i < cycleorpath.size()-1; i++) {
                        if(cycleorpath[i] == s ) {
                                idx = i;
                                break;
                        }
                }

                int entry = 0;
                if( idx == cycleorpath.size() - 1 ) {
                        if( cycleorpath[0] != s ) { 
                                entry = cycleorpath[0];
                        } else { // s -> .... -> s
                                entry = cycleorpath[1];
                        }
                } else {
                        entry = cycleorpath[idx + 1];
                }

                div_t block_tmp = div((int)entry,(int)config.k);
                PartitionID block = block_tmp.rem;
                unsigned number_of_vertices_to_move = block_tmp.quot + 1;

                for( unsigned i = 0; i < cycleorpath.size(); i++) {
                        div_t tmp = div((int)cycleorpath[i],(int)config.k);
                        if( (unsigned)tmp.rem == block && ((unsigned)(tmp.quot+1)>number_of_vertices_to_move)) {
                                number_of_vertices_to_move = tmp.quot + 1;
                        } 
                }


                if(boundary.getBlockWeight(block) == number_of_vertices_to_move) {
                        // then the block size would be zero afterwards
                        conflict_detected = true;
                }
        }


        return conflict_detected;

}
inline
bool advanced_models::handle_ultra_model_conflicts( PartitionConfig & config, 
                                                    graph_access & cycle_problem, 
                                                    complete_boundary & boundary,
                                                    std::vector<NodeID> & id_mapping, 
                                                    std::vector<bool> & feasable_edge, 
                                                    std::vector< NodeID > & cycle, 
                                                    NodeID & s,
                                                    augmented_Qgraph & aqg, 
                                                    bool remove_only_between_layers) {

        bool conflict_detected = cycleorpath_has_conflicts(config, boundary, cycle, s, aqg);
        if(conflict_detected) {
                NodeID blocked_edge = 0;
                if(remove_only_between_layers) {
                        std::vector<NodeID> eligible_idx;
                        for( unsigned i = 0; i < cycle.size()-2; i++) {
                                if( cycle[i] < s ) {
                                        div_t lhs_tmp = div((int)cycle[i],(int)config.k);
                                        div_t rhs_tmp = div((int)cycle[i+1],(int)config.k);
                                        if(lhs_tmp.quot != rhs_tmp.quot) {
                                                eligible_idx.push_back(i);
                                        }
                                }
                        }
                        blocked_edge  = random_functions::nextInt(0, eligible_idx.size() - 2);
                } else {
                        blocked_edge  = random_functions::nextInt(0, cycle.size() - 2);
                }

                forall_nodes(cycle_problem, node) {
                        forall_out_edges(cycle_problem, e, node) {
                                NodeID target = cycle_problem.getEdgeTarget(e);
                                if( cycle[blocked_edge] == node && cycle[blocked_edge+1] == target ) {
                                        feasable_edge[id_mapping[e]] = false;
                                        break;
                                }
                        } endfor
                } endfor

                conflicts++;
                return true;
        } 

        return false;
}

inline 
void advanced_models::perform_augmented_move( PartitionConfig & config, graph_access & G, 
                                              complete_boundary & boundary, 
                                              std::vector<NodeID> & cycleorpath, 
                                              NodeID & s, NodeID & t, 
                                              augmented_Qgraph & aqg) {

        for( unsigned i = 0; i < cycleorpath.size()-1; i++) {
                if(cycleorpath[i] == s || cycleorpath[i+1] == s || cycleorpath[i] == t || cycleorpath[i+1] == t) {
                        //removeing a path in the adv cycleorpath case
                        continue;
                }
                if( cycleorpath[i] == cycleorpath[i+1] + config.k || cycleorpath[i]+config.k == cycleorpath[i+1]) {

                        //augmented model case 
                        continue;
                }

                div_t lhs_tmp = div((int)cycleorpath[i],(int)config.k);
                div_t rhs_tmp = div((int)cycleorpath[i+1],(int)config.k);
                unsigned load_diff = lhs_tmp.quot+1; // this line should make the advanced expanded model feasable

                boundary_pair pair;
                pair.k   = config.k;
                pair.lhs = lhs_tmp.rem;
                pair.rhs = rhs_tmp.rem;

                std::vector<NodeID> vertices_of_move;
                std::vector<PartitionID> blocks_of_move;
                aqg.get_associated_vertices(pair, load_diff, vertices_of_move);
                aqg.get_associated_blocks(pair, load_diff, blocks_of_move);

                NodeWeight diff = 0;
                for( unsigned j = 0; j <  vertices_of_move.size(); j++) {
                        NodeID node = vertices_of_move[j];

                        PartitionID from = pair.lhs == blocks_of_move[j] ? pair.rhs : pair.lhs;
                        PartitionID to   = blocks_of_move[j];

                        G.setPartitionIndex(node, to);

                        boundary.postMovedBoundaryNodeUpdates(node, &pair, true, true);

                        NodeWeight this_nodes_weight = G.getNodeWeight(node);
                        boundary.setBlockNoNodes(from, boundary.getBlockNoNodes(from)-1);
                        boundary.setBlockNoNodes(to, boundary.getBlockNoNodes(to)+1);
                        boundary.setBlockWeight(from, boundary.getBlockWeight(from)-this_nodes_weight);
                        boundary.setBlockWeight(to, boundary.getBlockWeight(to)+this_nodes_weight);

                        diff += this_nodes_weight;
                }
        }
}
bool advanced_models::build_rebalance_model( PartitionConfig & config, 
                graph_access & G, 
                graph_access & G_bar,
                complete_boundary & boundary, 
                augmented_Qgraph & aqg, 
                std::vector<bool> & feasable_edge,
                unsigned & steps,
                graph_access & cycle_problem,
                NodeID & s,
                NodeID & t,
                std::vector<NodeID> & id_mapping) {

        //********************************************************************
        //build model 
        //********************************************************************
        unsigned max_vertex_weight_difference = aqg.get_max_vertex_weight_difference();
        NodeID number_of_nodes  = G_bar.number_of_nodes()*max_vertex_weight_difference + 2; 

        int max_diff_to_UB = 0;
        for( unsigned i = 0; i < config.k; i++) {
                int diff = config.upper_bound_partition - boundary.getBlockWeight(i);
                if(diff > max_diff_to_UB) {
                        max_diff_to_UB = diff;
                }
        }
        unsigned square = steps*(max_diff_to_UB+1);
        EdgeID number_of_edges  = G_bar.number_of_edges()*square + 4*number_of_nodes; 
        s = number_of_nodes - 2;
        t = number_of_nodes - 1;

        //initialize the deleted edges field
        if( feasable_edge.size() == 0 ) {
                feasable_edge.resize(number_of_edges);
                for( unsigned e = 0; e < feasable_edge.size(); e++) {
                        feasable_edge[e] = true;
                }
        }

        if( id_mapping.size() == 0 ) {
                id_mapping.resize( number_of_edges );
        }

        cycle_problem.start_construction(number_of_nodes, number_of_edges);
        unsigned edge_counter = 0;

        for( unsigned s_idx = 0; s_idx < max_vertex_weight_difference; s_idx++) {
                //create a new layer
                forall_nodes(G_bar, lhs) {
                        NodeID cur_node = cycle_problem.new_node();
                        forall_out_edges(G_bar, e, lhs) {
                                EdgeID rhs   = G_bar.getEdgeTarget(e);

                                //find the right edge in the augmented quotient graph
                                boundary_pair bp;
                                bp.k   = config.k;
                                bp.lhs = lhs;
                                bp.rhs = rhs;
                       
                                unsigned load_difference = s_idx + 1;
                                if( aqg.exists_vmovements_of_diff(bp, load_difference) ) {
                                        NodeID cur_target = s_idx*config.k+rhs;
                                        EdgeWeight cur_edge_weight = - aqg.get_gain_of_vmovements(bp, load_difference);

                                        if( feasable_edge[edge_counter] ) {
                                                EdgeID e_bar = cycle_problem.new_edge(cur_node, cur_target);
                                                cycle_problem.setEdgeWeight(e_bar, cur_edge_weight);
                                                id_mapping[e_bar] = edge_counter;
                                        }
                                        edge_counter++;

                                        //now the backward edges to the previous layers were possible
                                        NodeWeight cur_block_weight = boundary.getBlockWeight(rhs);
                                        if( s_idx != 0 ) {
                                                for( int j = s_idx, possible_overload = 1; j > 0 ; j--, possible_overload++) {
                                                        cur_target -= config.k;
                                                        if( cur_block_weight + possible_overload < config.upper_bound_partition) {
                                                                if( feasable_edge[edge_counter] ) {
                                                                        EdgeID e_bar = cycle_problem.new_edge(cur_node, cur_target);
                                                                        cycle_problem.setEdgeWeight(e_bar, cur_edge_weight);
                                                                        id_mapping[e_bar] = edge_counter;
                                                                }
                                                                edge_counter++;
                                                        } 
                                                }
                                        }

                                }
                        } endfor

                        //create a backward edge to t if it can take s_idx+1 vertices
                        if( boundary.getBlockWeight(lhs) + s_idx + 1 <= config.upper_bound_partition) {
                                if( feasable_edge[edge_counter] ) {
                                        EdgeID e_bar = cycle_problem.new_edge(cur_node, t);
                                        cycle_problem.setEdgeWeight(e_bar, 0);
                                }
                                edge_counter++;
                        }

                        if(s_idx != max_vertex_weight_difference-1) {
                                if( feasable_edge[edge_counter] ) {
                                        EdgeID e_bar = cycle_problem.new_edge(cur_node, cur_node+config.k);
                                        cycle_problem.setEdgeWeight(e_bar, 0);
                                        id_mapping[e_bar] = edge_counter;
                                }
                                edge_counter++;
                        }
                } endfor
        }

        s = cycle_problem.new_node();


        //now connect s to all vertices of the layer graph with weight zero
        for( unsigned s_idx = 0; s_idx < max_vertex_weight_difference; s_idx++) {
                forall_nodes(G_bar, node) {
                        if( boundary.getBlockWeight(node) > config.upper_bound_partition && feasable_edge[edge_counter]) {
                                EdgeID e = cycle_problem.new_edge(s, s_idx*config.k + node);
                                cycle_problem.setEdgeWeight(e, 0);
                                id_mapping[e] = edge_counter;
                        } 
                        edge_counter++;
                } endfor
        }

        t = cycle_problem.new_node();
        cycle_problem.finish_construction();

        return false;
}



//build ultra model ( removing deleted edges )
bool advanced_models::build_ultra_model( PartitionConfig & config, 
                                         graph_access & G, 
                                         graph_access & G_bar,
                                         complete_boundary & boundary, 
                                         augmented_Qgraph & aqg, 
                                         std::vector<bool> & feasable_edge,
                                         unsigned & steps,
                                         graph_access & cycle_problem,
                                         NodeID & s,
                                         std::vector<NodeID> & id_mapping) {

        //********************************************************************
        //build model 
        //********************************************************************
        unsigned max_vertex_weight_difference = aqg.get_max_vertex_weight_difference();
        NodeID number_of_nodes  = G_bar.number_of_nodes()*max_vertex_weight_difference+ 1; 

        int max_diff_to_UB = 0;
        for( unsigned i = 0; i < config.k; i++) {
                int diff = config.upper_bound_partition - boundary.getBlockWeight(i);
                if(diff > max_diff_to_UB) {
                        max_diff_to_UB = diff;
                }
        }
        unsigned square = steps*(max_diff_to_UB+1);
        EdgeID number_of_edges  = G_bar.number_of_edges()*square + 4*number_of_nodes; 
        s = number_of_nodes - 1;

        //initialize the deleted edges field
        if( feasable_edge.size() == 0 ) {
                feasable_edge.resize(number_of_edges);
                for( unsigned e = 0; e < feasable_edge.size(); e++) {
                        feasable_edge[e] = true;
                }
        }

        if( id_mapping.size() == 0 ) {
                id_mapping.resize( number_of_edges );
        }

        cycle_problem.start_construction(number_of_nodes, number_of_edges);
        unsigned edge_counter = 0;


        for( unsigned s_idx = 0; s_idx < max_vertex_weight_difference; s_idx++) {
                //create a new layer
                forall_nodes(G_bar, lhs) {
                        NodeID cur_node = cycle_problem.new_node();
                        forall_out_edges(G_bar, e, lhs) {
                                EdgeID rhs   = G_bar.getEdgeTarget(e);

                                //find the right edge in the augmented quotient graph
                                boundary_pair bp;
                                bp.k   = config.k;
                                bp.lhs = lhs;
                                bp.rhs = rhs;
                       
                                unsigned load_difference = s_idx + 1;
                                if( aqg.exists_vmovements_of_diff(bp, load_difference) ) {
                                        NodeID cur_target = s_idx*config.k+rhs;
                                        EdgeWeight cur_edge_weight = - aqg.get_gain_of_vmovements(bp, load_difference);

                                        if( feasable_edge[edge_counter] ) {
                                                EdgeID e_bar = cycle_problem.new_edge(cur_node, cur_target);
                                                cycle_problem.setEdgeWeight(e_bar, cur_edge_weight);
                                                id_mapping[e_bar] = edge_counter;
                                        }
                                        edge_counter++;

                                        // now the backward edges to the previous layers were possible
                                        NodeWeight cur_block_weight = boundary.getBlockWeight(rhs);
                                        if( s_idx != 0 ) {
                                                for( int j = s_idx, possible_overload = 1; j > 0 ; j--, possible_overload++) {
                                                        cur_target -= config.k;
                                                        if( cur_block_weight + possible_overload < config.upper_bound_partition) {
                                                                if( feasable_edge[edge_counter] ) {
                                                                        EdgeID e_bar = cycle_problem.new_edge(cur_node, cur_target);
                                                                        cycle_problem.setEdgeWeight(e_bar, cur_edge_weight);
                                                                        id_mapping[e_bar] = edge_counter;
                                                                }
                                                                edge_counter++;
                                                        } 
                                                }
                                        }

                                }
                        } endfor

                        if( boundary.getBlockWeight(lhs) + s_idx < config.upper_bound_partition) {
                                if( feasable_edge[edge_counter] ) {
                                        EdgeID e_bar = cycle_problem.new_edge(cur_node, s);
                                        cycle_problem.setEdgeWeight(e_bar, 0);
                                }
                                edge_counter++;
                        }
                        if(s_idx != max_vertex_weight_difference-1) {
                                if( feasable_edge[edge_counter] ) {
                                        EdgeID e_bar = cycle_problem.new_edge(cur_node, cur_node+config.k);
                                        id_mapping[e_bar] = edge_counter;
                                        cycle_problem.setEdgeWeight(e_bar, 0);
                                }
                                edge_counter++;
                        }
                } endfor
        }

        s = cycle_problem.new_node();

        //now connect s to all vertices of the layer graph with weight zero
        for( unsigned s_idx = 0; s_idx < max_vertex_weight_difference; s_idx++) {
                forall_nodes(G_bar, node) {
                        if( feasable_edge[edge_counter] ) {
                                EdgeID e = cycle_problem.new_edge(s, s_idx*config.k + node);
                                id_mapping[e] = edge_counter;
                                cycle_problem.setEdgeWeight(e, 0);
                        }
                        edge_counter++;
                } endfor
        }
        cycle_problem.finish_construction();

        return false;
}

#endif /* end of include guard: ADVANCED_MODELS_PR6SXN3G */
