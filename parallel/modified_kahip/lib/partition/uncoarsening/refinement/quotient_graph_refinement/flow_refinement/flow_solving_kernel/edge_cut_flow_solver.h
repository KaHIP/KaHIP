/******************************************************************************
 * edge_cut_flow_solver.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef EDGE_FLOW_SOLVER_4P49OMM
#define EDGE_FLOW_SOLVER_4P49OMM

#include "flow_solver.h"

class edge_cut_flow_solver : public flow_solver {
        public:
                edge_cut_flow_solver( );
                virtual ~edge_cut_flow_solver();

                EdgeWeight get_min_flow_max_cut(const PartitionConfig & config, 
                                                graph_access & G, 
                                                PartitionID & lhs, 
                                                PartitionID & rhs, 
                                                std::vector<NodeID> & lhs_boundary_stripe,
                                                std::vector<NodeID> & rhs_boundary_stripe,
                                                std::vector<NodeID> & new_to_old_ids,
                                                EdgeWeight & initial_cut,
                                                NodeWeight & rhs_part_weight,
                                                NodeWeight & rhs_stripe_weight,
                                                std::vector<NodeID> & new_rhs_nodes);               

                EdgeID regions_no_edges(graph_access & G,
                                        std::vector<NodeID> & lhs_boundary_stripe,
                                        std::vector<NodeID> & rhs_boundary_stripe,
                                        PartitionID & lhs, 
                                        PartitionID & rhs,
                                        std::vector<NodeID> & outer_lhs_boundary_nodes,
                                        std::vector<NodeID> & outer_rhs_boundary_nodes ); 


                //modified parse code from hi_pr
                EdgeWeight convert_ds(const PartitionConfig & config, 
                                      graph_access & G, 
                                      PartitionID & lhs, 
                                      PartitionID & rhs, 
                                      std::vector<NodeID> & lhs_boundary_stripe,
                                      std::vector<NodeID> & rhs_boundary_stripe,
                                      std::vector<NodeID> & new_to_old_ids,              
                                      long *n_ad, 
                                      long* m_ad, 
                                      node** nodes_ad, 
                                      arc** arcs_ad, 
                                      long ** cap_ad,
                                      node** source_ad, 
                                      node** sink_ad, 
                                      long* node_min_ad,
                                      EdgeID & no_edge_in_flow_graph); 

};




#endif /* end of include guard: FLOW_SOLVER_4P49OMM */
