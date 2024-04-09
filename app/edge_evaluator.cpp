/******************************************************************************
 * evaluator.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Daniel Seemaier <daniel.seemaier@kit.edu>
 *****************************************************************************/

#include <argtable3.h>
#include <regex.h>
#include <string.h>

#include "data_structure/graph_access.h"
#include "graph_io.h"
#include "parse_parameters.h"
#include "partition/partition_config.h"
#include "quality_metrics.h"

int main(int argn, char **argv) {
    PartitionConfig partition_config;
    std::string graph_filename;

    bool is_graph_weighted = false;
    bool suppress_output = false;
    bool recursive = false;

    int ret_code =
        parse_parameters(argn, argv, partition_config, graph_filename,
                         is_graph_weighted, suppress_output, recursive);

    if (ret_code) {
        return 0;
    }

    graph_access G;
    graph_io::readGraphWeighted(G, graph_filename);
    G.set_partition_count(partition_config.k);

    std::vector<EdgeID> edge_partition(G.number_of_edges());

    std::vector<PartitionID> input_partition;
    if (partition_config.input_partition != "") {
        std::cout << "reading input partition" << std::endl;
        graph_io::readVector(edge_partition, partition_config.input_partition);
    } else {
        std::cout << "Please specify an input partition using the "
                     "--input_partition flag."
                  << std::endl;
        exit(0);
    }

    std::cout << "graph has " << G.number_of_nodes() << " nodes and "
              << G.number_of_edges() << " edges" << std::endl;

    for (NodeID u = 0; u < G.number_of_nodes(); ++u) {
        for (EdgeID e = G.get_first_edge(u); e < G.get_first_invalid_edge(u);
             ++e) {
            NodeID v = G.getEdgeTarget(e);
            for (EdgeID e_prime = G.get_first_edge(v);
                 e_prime < G.get_first_invalid_edge(v); ++e_prime) {
                if (G.getEdgeTarget(e_prime) == u) {
                    PartitionID part = edge_partition[e];
                    PartitionID part_prime = edge_partition[e_prime];
                    if (part != part_prime) {
                        std::cerr << "edge (" << u << ", " << v
                                  << ") has different partitions: " << part
                                  << " and " << part_prime << std::endl;
                        std::exit(1);
                    }
                }
            }
        }
    }

    unsigned cost = 0;
    std::vector<bool> counted(G.get_partition_count());
    for (NodeID u = 0; u < G.number_of_nodes(); ++u) {
        if (G.getNodeDegree(u) == 0) continue;

        for (EdgeID e = G.get_first_edge(u); e < G.get_first_invalid_edge(u);
             ++e) {
            PartitionID part = edge_partition[e];
            if (!counted[part]) {
                counted[part] = true;
                ++cost;
            }
        }

        --cost;
        std::fill(counted.begin(), counted.end(), false);
    }

    std::cout << "vertex cut: " << cost << std::endl;

    quality_metrics qm;
    std::cout << "balance: " << qm.edge_balance(G, edge_partition) << std::endl;
}

