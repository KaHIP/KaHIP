/******************************************************************************
 * dspac.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Graph Partitioning
 ******************************************************************************
 * Copyright (C) 2018 Christian Schulz
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
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
