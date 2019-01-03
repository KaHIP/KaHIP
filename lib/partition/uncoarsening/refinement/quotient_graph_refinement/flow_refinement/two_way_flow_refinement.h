/******************************************************************************
 * two_way_flow_refinement.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef TWO_WAY_FLOW_REFINEMENT_BVTL6G49
#define TWO_WAY_FLOW_REFINEMENT_BVTL6G49

#include "data_structure/graph_access.h"
#include "data_structure/priority_queues/priority_queue_interface.h"
#include "partition_config.h"
#include "uncoarsening/refinement/quotient_graph_refinement/complete_boundary.h"
#include "uncoarsening/refinement/quotient_graph_refinement/partial_boundary.h"
#include "uncoarsening/refinement/quotient_graph_refinement/two_way_refinement.h"

class two_way_flow_refinement : public two_way_refinement {
        public:
                two_way_flow_refinement( );
                virtual ~two_way_flow_refinement();

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
        private:
                EdgeWeight iterativ_flow_iteration(PartitionConfig & config, 
                                                   graph_access & G,
                                                   complete_boundary & boundary, 
                                                   std::vector<NodeID> & lhs_pq_start_nodes, 
                                                   std::vector<NodeID> & rhs_pq_start_nodes,
                                                   boundary_pair * refinement_pair,        
                                                   NodeWeight & lhs_part_weight,
                                                   NodeWeight & rhs_part_weight,
                                                   EdgeWeight & cut,
                                                   bool & something_changed); 

                void apply_partition_and_update_boundary(const PartitionConfig & config, 
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
                                                         std::vector<NodeID> & new_rhs_nodes); 


};


#endif /* end of include guard: TWO_WAY_FLOW_REFINEMENT_BVTL6G49 */
