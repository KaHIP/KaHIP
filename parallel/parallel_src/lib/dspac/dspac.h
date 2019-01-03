/******************************************************************************
 * dspac.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Author: Daniel Seemaier <daniel.seemaier@student.kit.edu>
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef KAHIP_DSPAC_H
#define KAHIP_DSPAC_H

#include <vector>
#include "data_structure/parallel_graph_access.h"
#include "definitions.h"

class dspac {
public:
    dspac(parallel_graph_access &graph, MPI_Comm comm, EdgeWeight infinity);
    void construct(parallel_graph_access &split_graph);
    std::vector<PartitionID> project_partition(parallel_graph_access &split_graph);
    EdgeWeight calculate_vertex_cut(PartitionID k, const std::vector<PartitionID> &edge_partition);
    void fix_cut_dominant_edges(parallel_graph_access &split_graph);

private:
    bool assert_adjacency_lists_sorted();
    bool assert_sanity_checks(parallel_graph_access &split_graph);
    bool assert_edge_range_array_ok(const std::vector<NodeID> &edge_range_array);
    bool assert_node_range_array_ok(const std::vector<NodeID> &node_range_array);

    void internal_construct(parallel_graph_access &split_graph);

    MPI_Comm m_comm;
    EdgeWeight m_infinity;
    parallel_graph_access &m_input_graph;
};

#endif // KAHIP_DSPAC_H
