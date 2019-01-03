/******************************************************************************
 * vertex_separator_algorithm.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef VERTEX_SEPARTATOR_ALGORITHM_XUDNZZM8
#define VERTEX_SEPARTATOR_ALGORITHM_XUDNZZM8

#include <unordered_map>

#include "data_structure/graph_access.h"
#include "partition_config.h"
#include "uncoarsening/refinement/quotient_graph_refinement/complete_boundary.h"

class vertex_separator_algorithm {
        public:
                vertex_separator_algorithm();
                virtual ~vertex_separator_algorithm();

                void compute_vertex_separator(const PartitionConfig & config, 
                                              graph_access & G, 
                                              complete_boundary & boundary, 
                                              std::vector<NodeID> & overall_separator);

                void compute_vertex_separator_simple(const PartitionConfig & config, 
                                                     graph_access & G, 
                                                     complete_boundary & boundary, 
                                                     std::vector<NodeID> & overall_separator);

                void compute_vertex_separator(const PartitionConfig & config, 
                                              graph_access & G, 
                                              complete_boundary & boundary);

                //ASSERTIONS
                bool is_vertex_separator(graph_access & G, std::unordered_map<NodeID, bool> & separator);

};


#endif /* end of include guard: VERTEX_SEPARTATOR_ALGORITHM_XUDNZZM8 */
