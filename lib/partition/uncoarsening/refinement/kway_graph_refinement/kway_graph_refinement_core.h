/******************************************************************************
 * kway_graph_refinement_core.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef KWAY_GRAPH_REFINEMENT_CORE_PVGY97EW
#define KWAY_GRAPH_REFINEMENT_CORE_PVGY97EW

#include <unordered_map>
#include <vector>

#include "data_structure/priority_queues/priority_queue_interface.h"
#include "definitions.h"
#include "kway_graph_refinement_commons.h"
#include "tools/random_functions.h"
#include "uncoarsening/refinement/quotient_graph_refinement/2way_fm_refinement/vertex_moved_hashtable.h"
#include "uncoarsening/refinement/refinement.h"

class kway_graph_refinement_core  {
        public:
                kway_graph_refinement_core( );
                virtual ~kway_graph_refinement_core();

                EdgeWeight single_kway_refinement_round(PartitionConfig & config, 
                                                        graph_access & G, 
                                                        complete_boundary & boundary, 
                                                        boundary_starting_nodes & start_nodes, 
                                                        int step_limit, 
                                                        vertex_moved_hashtable & moved_idx );

                EdgeWeight single_kway_refinement_round(PartitionConfig & config, 
                                                        graph_access & G, 
                                                        complete_boundary & boundary, 
                                                        boundary_starting_nodes & start_nodes, 
                                                        int step_limit, 
                                                        vertex_moved_hashtable & moved_idx,
                                                        std::unordered_map<PartitionID, PartitionID> & touched_blocks); 


         private:
               EdgeWeight single_kway_refinement_round_internal(PartitionConfig & config, 
                                                                graph_access & G, 
                                                                complete_boundary & boundary, 
                                                                boundary_starting_nodes & start_nodes, 
                                                                int step_limit,
                                                                vertex_moved_hashtable & moved_idx,
                                                                bool compute_touched_partitions,
                                                                std::unordered_map<PartitionID, PartitionID> &  touched_blocks); 


                void init_queue_with_boundary(const PartitionConfig & config, 
                                              graph_access & G, 
                                              std::vector<NodeID> & bnd_nodes, 
                                              refinement_pq * queue, 
                                              vertex_moved_hashtable & moved_idx);

                inline bool move_node(PartitionConfig & config, 
                                      graph_access & G, 
                                      NodeID & node, 
                                      vertex_moved_hashtable & moved_idx, 
                                      refinement_pq * queue, 
                                      complete_boundary & boundary);

                inline void move_node_back(PartitionConfig & config, 
                                           graph_access & G, 
                                           NodeID & node,
                                           PartitionID & to, 
                                           vertex_moved_hashtable & moved_idx, 
                                           refinement_pq * queue, 
                                           complete_boundary & boundary);

                void initialize_partition_moves_array(PartitionConfig & config, 
                                                      complete_boundary & boundary, 
                                                      std::vector<bool> & partition_move_valid); 
                
                kway_graph_refinement_commons* commons;
};

inline bool kway_graph_refinement_core::move_node(PartitionConfig & config, 
                graph_access & G, 
                NodeID & node, 
                vertex_moved_hashtable & moved_idx, 
                refinement_pq * queue, 
                complete_boundary & boundary) {

        ASSERT_TRUE(boundary.assert_bnodes_in_boundaries());
        ASSERT_TRUE(boundary.assert_boundaries_are_bnodes());

        PartitionID from = G.getPartitionIndex(node);
        PartitionID to;
        EdgeWeight node_ext_deg;
        commons->compute_gain(G, node, to, node_ext_deg);

        NodeWeight this_nodes_weight = G.getNodeWeight(node);
        if(boundary.getBlockWeight(to) + this_nodes_weight >= config.upper_bound_partition) 
                return false;

        if(boundary.getBlockNoNodes(from) - 1 == 0) // assure that no block gets accidentally empty
                return false;

        G.setPartitionIndex(node, to);        

        boundary_pair pair;
        pair.k = config.k;
        pair.lhs = from;
        pair.rhs = to;

        boundary.postMovedBoundaryNodeUpdates(node, &pair, true, true);

        boundary.setBlockNoNodes(from, boundary.getBlockNoNodes(from)-1);
        boundary.setBlockNoNodes(to,   boundary.getBlockNoNodes(to)+1);
        boundary.setBlockWeight( from, boundary.getBlockWeight(from)-this_nodes_weight);
        boundary.setBlockWeight( to,   boundary.getBlockWeight(to)+this_nodes_weight);

        ASSERT_TRUE(boundary.assert_bnodes_in_boundaries());
        ASSERT_TRUE(boundary.assert_boundaries_are_bnodes());

        //update gain of neighbors / the boundaries have allready been updated
        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                PartitionID targets_max_gainer;
                EdgeWeight ext_degree; // the local external degree
                Gain gain = commons->compute_gain(G, target, targets_max_gainer, ext_degree);

                if(queue->contains(target)) {
                        assert(moved_idx.find(target) != moved_idx.end());
                        if(ext_degree > 0) {
                                queue->changeKey(target, gain);
                        } else {
                                queue->deleteNode(target);
                        }
                } else {
                        if(ext_degree > 0) {
                                if(moved_idx.find(target) == moved_idx.end()) {
                                        queue->insert(target, gain);
                                        moved_idx[target].index = NOT_MOVED;
                                } 
                        } 
                }
        } endfor

        return true;
}

#endif //[> end of include guard: KWAY_GRAPH_REFINEMENT_PVGY97EW <]

