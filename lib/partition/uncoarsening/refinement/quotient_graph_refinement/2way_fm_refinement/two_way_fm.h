/******************************************************************************
 * two_way_fm.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef TWOWAY_FM_YLYN82Y1
#define TWOWAY_FM_YLYN82Y1

#include <vector>

#include "data_structure/graph_access.h"
#include "data_structure/priority_queues/priority_queue_interface.h"
#include "definitions.h"
#include "partition_config.h"
#include "uncoarsening/refinement/quotient_graph_refinement/complete_boundary.h"
#include "uncoarsening/refinement/quotient_graph_refinement/partial_boundary.h"
#include "uncoarsening/refinement/quotient_graph_refinement/two_way_refinement.h"
#include "vertex_moved_hashtable.h"


class two_way_fm : public two_way_refinement {
        public:
                two_way_fm( );
                virtual ~two_way_fm();
                EdgeWeight perform_refinement(PartitionConfig & config, 
                                graph_access & G, 
                                complete_boundary & boundary, 
                                std::vector<NodeID> & lhs_pq_start_nodes, 
                                std::vector<NodeID> & rhs_pq_start_nodes,
                                boundary_pair * refinement_pair,        
                                NodeWeight & lhs_part_weight,
                                NodeWeight & rhs_part_weight,
                                EdgeWeight & cut,
                                bool & something_changed);

                inline bool int_ext_degree(graph_access & G, 
                                const NodeID & node,
                                const PartitionID lhs,
                                const PartitionID rhs,
                                EdgeWeight & int_degree,
                                EdgeWeight & ext_degree);


        private:
                void init_queue_with_boundary(const PartitionConfig & config,
                                graph_access & G,
                                std::vector<NodeID> &bnd_nodes,
                                refinement_pq * queue,                     
                                PartitionID partition_of_boundary, 
                                PartitionID other); 


                void move_node(const PartitionConfig & config, 
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
                               complete_boundary & boundary);

                void move_node_back(const PartitionConfig & config, 
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
                                    complete_boundary & boundary); 


                ///////////////////////////////////////////////////////////////////////////
                //Assertions
                ///////////////////////////////////////////////////////////////////////////
#ifndef NDEBUG
                //assert that every node in the lhs boundary has external degree > 0
                bool assert_only_boundary_nodes(graph_access & G, 
                                                PartialBoundary & lhs_boundary, 
                                                PartitionID lhs, 
                                                PartitionID rhs);

                //assert that every node with ext degree > 0 is lhs boundary 
                bool assert_every_boundary_nodes(graph_access & G, 
                                                 PartialBoundary & lhs_boundary, 
                                                 PartitionID lhs, 
                                                 PartitionID rhs);

                //check all of the possible compinations of the two assertions above
                bool assert_directed_boundary_condition(graph_access & G, 
                                                        complete_boundary & boundary, 
                                                        PartitionID lhs, 
                                                        PartitionID rhs);
#endif

};

inline bool two_way_fm::int_ext_degree( graph_access & G, 
                const NodeID & node,
                const PartitionID lhs,
                const PartitionID rhs,
                EdgeWeight & int_degree,
                EdgeWeight & ext_degree) {


        ASSERT_EQ(lhs, G.getPartitionIndex(node));

        int_degree               = 0;
        ext_degree               = 0;
        bool update_is_difficult = false;

        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                PartitionID targets_partition = G.getPartitionIndex(target);

                if(targets_partition == lhs) {
                        int_degree += G.getEdgeWeight(e); 
                } else if(targets_partition == rhs) {
                        ext_degree += G.getEdgeWeight(e);
                }

                if(targets_partition != lhs && targets_partition != rhs) {
                        update_is_difficult = true;
                } 
        } endfor

        return update_is_difficult;
}


#endif /* end of include guard: TWO_WAY_FM_YLYN82Y1 */
