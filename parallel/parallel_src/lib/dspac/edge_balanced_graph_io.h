/******************************************************************************
 * edge_balanced_graph_io.h
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

#ifndef KAHIP_EDGEBALANCED_GRAPH_IO_H
#define KAHIP_EDGEBALANCED_GRAPH_IO_H

#include <fstream>
#include <iostream>

#include "definitions.h"
#include "data_structure/parallel_graph_access.h"
#include "partition_config.h"

class edge_balanced_graph_io {
public:
    static void read_binary_graph_edge_balanced(parallel_graph_access &G, const std::string &filename,
            const PPartitionConfig &config, int rank, int size);

    static void read_binary_graph_edge_balanced(parallel_graph_access &G, const std::string &filename,
            const PPartitionConfig &config);
};

#endif // KAHIP_EDGEBALANCED_GRAPH_IO_H