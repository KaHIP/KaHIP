/******************************************************************************
 * quality_metrics.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/


#ifndef QUALITY_METRICS_10HC2I5M
#define QUALITY_METRICS_10HC2I5M

#include "data_structure/graph_access.h"
#include "data_structure/matrix/matrix.h"
#include "partition_config.h"

class quality_metrics {
public:
        quality_metrics();
        virtual ~quality_metrics ();

        EdgeWeight edge_cut(graph_access & G);
        EdgeWeight edge_cut(graph_access & G, int * partition_map); 
        EdgeWeight edge_cut(graph_access & G, PartitionID lhs, PartitionID rhs);
        EdgeWeight max_communication_volume(graph_access & G);
        EdgeWeight min_communication_volume(graph_access & G);
        EdgeWeight max_communication_volume(graph_access & G, int * partition_map);
        EdgeWeight total_communication_volume(graph_access & G); 
        EdgeWeight objective(const PartitionConfig & config, graph_access & G, int * partition_map);
        EdgeWeight edge_cut_connected(graph_access & G, int * partition_map);
        int boundary_nodes(graph_access & G);
        NodeWeight separator_weight(graph_access& G);
        double balance(graph_access & G);
        double balance_edges(graph_access & G);
        double balance_separator(graph_access & G);
        double edge_balance(graph_access &G, const std::vector<PartitionID> &edge_partition);

        NodeWeight total_qap(graph_access & C, matrix & D, std::vector< NodeID > & rank_assign);
        NodeWeight total_qap(matrix & C, matrix & D, std::vector< NodeID > & rank_assign);
};

#endif /* end of include guard: QUALITY_METRICS_10HC2I5M */
