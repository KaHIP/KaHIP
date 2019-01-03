/******************************************************************************
 * spac.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Author: Seemaier, Daniel Max Manfred <daniel.seemaier@student.kit.edu>
 *****************************************************************************/

#ifndef KAHIP_SPAC_H
#define KAHIP_SPAC_H

#include <vector>

#include "data_structure/graph_access.h"
#include "definitions.h"

class spac {
public:
    spac(graph_access &input_graph, EdgeWeight infinity);

    graph_access &construct_split_graph();

    void fix_cut_dominant_edges();

    std::vector<PartitionID> project_partition();

    unsigned calculate_vertex_cut(const std::vector<PartitionID> &edge_partition);

private:
    void find_reverse_edges();

    graph_access &m_input_graph;
    graph_access m_split_graph;
    EdgeWeight m_infinity;
    std::vector<EdgeID> m_reverse_edge;
};

#endif // KAHIP_SPAC_H
