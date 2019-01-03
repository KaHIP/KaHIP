/******************************************************************************
 * two_way_flow_refinement.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <sstream>

#include "boundary_bfs.h"
#include "flow_solving_kernel/cut_flow_problem_solver.h"
#include "quality_metrics.h"
#include "two_way_flow_refinement.h"

two_way_flow_refinement::two_way_flow_refinement() {

}

two_way_flow_refinement::~two_way_flow_refinement() {

}

EdgeWeight two_way_flow_refinement::perform_refinement(PartitionConfig & config, 
                                                       graph_access & G, 
                                                       complete_boundary & boundary, 
                                                       std::vector<NodeID> & lhs_pq_start_nodes, 
                                                       std::vector<NodeID> & rhs_pq_start_nodes,
                                                       boundary_pair * refinement_pair,        
                                                       NodeWeight & lhs_part_weight,
                                                       NodeWeight & rhs_part_weight,
                                                       EdgeWeight & cut,
                                                       bool & something_changed) {

        EdgeWeight retval  =  iterativ_flow_iteration(config, G, boundary, lhs_pq_start_nodes, rhs_pq_start_nodes, 
                                                      refinement_pair, lhs_part_weight, rhs_part_weight, cut, something_changed);

        if(retval > 0) {
                something_changed = true;
        }

        return retval;
}


EdgeWeight two_way_flow_refinement::iterativ_flow_iteration(PartitionConfig & config, 
                                                            graph_access & G,
                                                            complete_boundary & boundary, 
                                                            std::vector<NodeID> & lhs_pq_start_nodes, 
                                                            std::vector<NodeID> & rhs_pq_start_nodes,
                                                            boundary_pair * refinement_pair,        
                                                            NodeWeight & lhs_part_weight,
                                                            NodeWeight & rhs_part_weight,
                                                            EdgeWeight & cut,
                                                            bool & something_changed) {

        if(lhs_pq_start_nodes.size() == 0 or rhs_pq_start_nodes.size() == 0) return 0; // nothing to refine
        ASSERT_TRUE(lhs_part_weight < config.upper_bound_partition && rhs_part_weight < config.upper_bound_partition);

        PartitionID lhs = refinement_pair->lhs;
        PartitionID rhs = refinement_pair->rhs;
        boundary_bfs bfs_region_searcher;               

        double region_factor    = config.flow_region_factor;
        unsigned max_iterations = config.max_flow_iterations;
        unsigned iteration      = 0;

        std::vector<NodeID> lhs_nodes;
        std::vector<NodeID> rhs_nodes;

        EdgeWeight cur_improvement = 1;
        EdgeWeight best_cut = cut; 
        bool sumoverweight = lhs_part_weight + rhs_part_weight > 2*config.upper_bound_partition;
        if(sumoverweight) {
                return 0; 
        }

        NodeWeight average_partition_weight = ceil(config.work_load / config.k);
        while(cur_improvement > 0 && iteration < max_iterations) {
                NodeWeight upper_bound_no_lhs = (NodeWeight)std::max((100.0+region_factor*config.imbalance)/100.0*(average_partition_weight) - rhs_part_weight,0.0);
                NodeWeight upper_bound_no_rhs = (NodeWeight)std::max((100.0+region_factor*config.imbalance)/100.0*(average_partition_weight) - lhs_part_weight,0.0);

                upper_bound_no_lhs = std::min( lhs_part_weight-1, upper_bound_no_lhs);
                upper_bound_no_rhs = std::min( rhs_part_weight-1, upper_bound_no_rhs);

                std::vector<NodeID> lhs_boundary_stripe;
                NodeWeight lhs_stripe_weight = 0;
                if(!bfs_region_searcher.boundary_bfs_search(G, lhs_pq_start_nodes, lhs, 
                                        upper_bound_no_lhs, lhs_boundary_stripe, 
                                        lhs_stripe_weight, true)) {

                        EdgeWeight improvement = cut-best_cut;
                        cut = best_cut;
                        return improvement;
                }


                std::vector<NodeID> rhs_boundary_stripe;
                NodeWeight rhs_stripe_weight = 0;
                if(!bfs_region_searcher.boundary_bfs_search(G, rhs_pq_start_nodes, rhs, 
                                        upper_bound_no_rhs, rhs_boundary_stripe, 
                                        rhs_stripe_weight, true)) { 

                        EdgeWeight improvement = cut-best_cut;
                        cut = best_cut;
                        return improvement;
                }

                std::vector<NodeID> new_rhs_nodes;
                std::vector<NodeID> new_to_old_ids;

                cut_flow_problem_solver fsolve;
                EdgeWeight new_cut = fsolve.get_min_flow_max_cut(config, G, 
                                                                lhs, rhs, 
                                                                lhs_boundary_stripe, rhs_boundary_stripe, 
                                                                new_to_old_ids, best_cut, 
                                                                rhs_part_weight,
                                                                rhs_stripe_weight,
                                                                new_rhs_nodes);

                NodeWeight new_lhs_part_weight   = 0;
                NodeWeight new_rhs_part_weight   = 0;
                NodeWeight new_lhs_stripe_weight = 0;
                NodeWeight new_rhs_stripe_weight = 0;
                NodeID no_nodes_flow_graph = lhs_boundary_stripe.size() + rhs_boundary_stripe.size();

                for(unsigned i = 0; i < new_rhs_nodes.size(); i++) {
                        NodeID new_rhs_node = new_rhs_nodes[i];
                        if(new_rhs_node < no_nodes_flow_graph) { // not target and source
                                NodeID old_node_id = new_to_old_ids[new_rhs_node];
                                new_rhs_stripe_weight += G.getNodeWeight(old_node_id);
                                G.setPartitionIndex(old_node_id, BOUNDARY_STRIPE_NODE-1);
                        }                
                }

                for(unsigned i = 0; i < lhs_boundary_stripe.size(); i++) {
                        if(  G.getPartitionIndex(lhs_boundary_stripe[i]) == BOUNDARY_STRIPE_NODE) { 
                                new_lhs_stripe_weight += G.getNodeWeight(lhs_boundary_stripe[i]);
                        }
                }

                for(unsigned i = 0; i < rhs_boundary_stripe.size(); i++) {
                        if(  G.getPartitionIndex(rhs_boundary_stripe[i]) == BOUNDARY_STRIPE_NODE) { 
                                new_lhs_stripe_weight += G.getNodeWeight(rhs_boundary_stripe[i]);
                        }
                }
                new_lhs_part_weight = boundary.getBlockWeight(lhs) + ((int)new_lhs_stripe_weight-(int)lhs_stripe_weight) ;
                new_rhs_part_weight = boundary.getBlockWeight(rhs) + ((int)new_rhs_stripe_weight-(int)rhs_stripe_weight) ;

                bool partition_is_feasable = false;
                if(config.most_balanced_minimum_cuts) {
                       partition_is_feasable = new_lhs_part_weight < config.upper_bound_partition 
                                             && new_rhs_part_weight < config.upper_bound_partition 
                                             && (new_cut < best_cut || abs((int)new_lhs_part_weight - (int) new_rhs_part_weight) < abs((int)lhs_part_weight - (int)rhs_part_weight));
                } else {
                       partition_is_feasable = new_lhs_part_weight < config.upper_bound_partition 
                                            && new_rhs_part_weight < config.upper_bound_partition && new_cut < best_cut;
                }
                
                if(partition_is_feasable) {
                        // then this partition can be accepted
                        apply_partition_and_update_boundary( config, G, refinement_pair, 
                                                             lhs, rhs, boundary,
                                                             lhs_boundary_stripe, rhs_boundary_stripe,
                                                             lhs_stripe_weight, 
                                                             rhs_stripe_weight, 
                                                             new_to_old_ids,
                                                             new_rhs_nodes); 

                        boundary.setEdgeCut(refinement_pair, new_cut);

                        lhs_part_weight = boundary.getBlockWeight(lhs);
                        rhs_part_weight = boundary.getBlockWeight(rhs);
                        ASSERT_TRUE(lhs_part_weight < config.upper_bound_partition && rhs_part_weight < config.upper_bound_partition);

                       
                        cur_improvement = best_cut - new_cut;
                        best_cut        = new_cut;

                        if(2*region_factor < config.flow_region_factor) {
                                region_factor *= 2; 
                        } else {
                                region_factor = config.flow_region_factor;
                        }
                        if(region_factor == config.flow_region_factor) {
                                //in that case we are finished
                                break;
                        }
                        if( iteration+1 < max_iterations ) {
                                // update the start nodes for the bfs
                                lhs_pq_start_nodes.clear();
                                boundary.setup_start_nodes(G, lhs, *refinement_pair, lhs_pq_start_nodes); 

                                rhs_pq_start_nodes.clear(); 
                                boundary.setup_start_nodes(G, rhs, *refinement_pair, rhs_pq_start_nodes);
                        }
                        
                        
                } else {
                        //undo changes
                        for( unsigned i = 0; i < lhs_boundary_stripe.size(); i++) {
                                G.setPartitionIndex(lhs_boundary_stripe[i], lhs); 
                        }
                        for( unsigned i = 0; i < rhs_boundary_stripe.size(); i++) {
                                G.setPartitionIndex(rhs_boundary_stripe[i], rhs); 
                        }

                        //smaller the region_factor
                        region_factor = std::max(region_factor/2,1.0);
                        if(new_cut == best_cut) {
                                break; 
                        }
                }
                iteration++;
        }


        ASSERT_TRUE(lhs_part_weight < config.upper_bound_partition && rhs_part_weight < config.upper_bound_partition);
        EdgeWeight improvement = cut-best_cut;
        cut = best_cut;
        return improvement;
}


void two_way_flow_refinement::apply_partition_and_update_boundary( const PartitionConfig & config, 
                                                                   graph_access & G, 
                                                                   boundary_pair * refinement_pair,
                                                                   PartitionID & lhs, 
                                                                   PartitionID & rhs,
                                                                   complete_boundary & boundary, 
                                                                   std::vector<NodeID> & lhs_boundary_stripe,
                                                                   std::vector<NodeID> & rhs_boundary_stripe,
                                                                   NodeWeight & lhs_stripe_weight, 
                                                                   NodeWeight & rhs_stripe_weight, 
                                                                   std::vector<NodeID> & new_to_old_ids,
                                                                   std::vector<NodeID> & new_rhs_nodes) {


        NodeID no_nodes_flow_graph = lhs_boundary_stripe.size() + rhs_boundary_stripe.size();
        NodeWeight new_lhs_stripe_weight = 0;
        NodeWeight new_rhs_stripe_weight = 0;

        NodeWeight new_rhs_stripe_no_nodes = 0; 
        NodeWeight new_lhs_stripe_no_nodes = 0;

        for(unsigned i = 0; i < new_rhs_nodes.size(); i++) {
                NodeID new_rhs_node = new_rhs_nodes[i];
                if(new_rhs_node < no_nodes_flow_graph) { // not target and source
                        NodeID old_node_id = new_to_old_ids[new_rhs_node];
                        G.setPartitionIndex(old_node_id, rhs);
                        new_rhs_stripe_weight += G.getNodeWeight(old_node_id);
                        new_rhs_stripe_no_nodes++;
                }                
        }

        for(unsigned i = 0; i < lhs_boundary_stripe.size(); i++) {
                if(  G.getPartitionIndex(lhs_boundary_stripe[i]) == BOUNDARY_STRIPE_NODE) { 
                        G.setPartitionIndex(lhs_boundary_stripe[i], lhs);	
                        new_lhs_stripe_weight += G.getNodeWeight(lhs_boundary_stripe[i]);
                        new_lhs_stripe_no_nodes++;
                }
        }

        for(unsigned i = 0; i < rhs_boundary_stripe.size(); i++) {
                if(  G.getPartitionIndex(rhs_boundary_stripe[i]) == BOUNDARY_STRIPE_NODE) { 
                        G.setPartitionIndex(rhs_boundary_stripe[i], lhs);	
                        new_lhs_stripe_weight += G.getNodeWeight(rhs_boundary_stripe[i]);
                        new_lhs_stripe_no_nodes++;
                }
        }


        // ********** fix the boundary data structure *****************  
        boundary.setBlockWeight(lhs, boundary.getBlockWeight(lhs) + ((int)new_lhs_stripe_weight-(int)lhs_stripe_weight) );
        boundary.setBlockWeight(rhs, boundary.getBlockWeight(rhs) + ((int)new_rhs_stripe_weight-(int)rhs_stripe_weight) );

        boundary.setBlockNoNodes(lhs, boundary.getBlockNoNodes(lhs) + ((int)new_lhs_stripe_no_nodes -(int)lhs_boundary_stripe.size()) );
        boundary.setBlockNoNodes(rhs, boundary.getBlockNoNodes(rhs) + ((int)new_rhs_stripe_no_nodes -(int)rhs_boundary_stripe.size()) );


        //this can be improved by only calling this method on the nodes that changed the partition
        for(unsigned i = 0; i < lhs_boundary_stripe.size(); i++) {
                boundary.postMovedBoundaryNodeUpdates(lhs_boundary_stripe[i], refinement_pair, false, true); 
        }

        for(unsigned i = 0; i < rhs_boundary_stripe.size(); i++) {
                boundary.postMovedBoundaryNodeUpdates(rhs_boundary_stripe[i], refinement_pair, false, true); 
        }

} 


