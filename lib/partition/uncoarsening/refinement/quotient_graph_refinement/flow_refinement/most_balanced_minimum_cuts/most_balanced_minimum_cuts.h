/******************************************************************************
 * most_balanced_minimum_cuts.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef MOST_BALANCED_MINIMUM_CUTS_SBD5CS
#define MOST_BALANCED_MINIMUM_CUTS_SBD5CS

#include "data_structure/graph_access.h"
#include "partition_config.h"

class most_balanced_minimum_cuts {
        public:
                most_balanced_minimum_cuts();
                virtual ~most_balanced_minimum_cuts();

                void compute_good_balanced_min_cut( graph_access & residualGraph,
                                                    const PartitionConfig & config,
                                                    NodeWeight & perfect_rhs_weight, 
                                                    std::vector< NodeID > & new_rhs_node ); 

                // the same on a graph given as arrays: the edges of node v are
                // targets[first[v] .. first[v+1]), its weight is weights[v]; the last two
                // nodes are the source and the sink
                void compute_good_balanced_min_cut( const std::vector<EdgeID> & first,
                                                    const std::vector<NodeID> & targets,
                                                    const std::vector<NodeWeight> & weights,
                                                    const PartitionConfig & config,
                                                    NodeWeight & perfect_rhs_weight, 
                                                    std::vector< NodeID > & new_rhs_node ); 

        private:
                void compute_new_rhs( graph_access & scc_graph, 
                                      const PartitionConfig & config,
                                      std::vector< NodeWeight > & comp_weights,
                                      int comp_of_s,
                                      int comp_of_t,
                                      NodeWeight optimal_rhs_weight,
                                      std::vector<int> & comp_for_rhs); 
};


#endif /* end of include guard: MOST_BALANCED_MINIMUM_CUTS_SBD5CS */
