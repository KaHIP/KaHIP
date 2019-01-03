/******************************************************************************
 * cut_flow_problem_solver.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef CUT_FLOW_PROBLEM_SOLVER_4P49OMM
#define CUT_FLOW_PROBLEM_SOLVER_4P49OMM

#include "partition_config.h"
#include "data_structure/flow_graph.h"

class cut_flow_problem_solver  {
        public:
                cut_flow_problem_solver( );
                virtual ~cut_flow_problem_solver();

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
                                      flow_graph & rG); 

};




#endif /* end of include guard: FLOW_SOLVER_4P49OMM */
