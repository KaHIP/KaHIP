/******************************************************************************
 * edge_balanced_graph_io.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
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
