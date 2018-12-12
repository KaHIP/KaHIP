/******************************************************************************
 * spac.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2018
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
