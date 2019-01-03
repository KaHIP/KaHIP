/******************************************************************************
 * spac.cpp
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

#include <argtable2.h>
#include <iostream>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h>

#include "balance_configuration.h"
#include "data_structure/graph_access.h"
#include "parse_spac_parameters.h"
#include "partition/partition_config.h"
#include "partition/graph_partitioner.h"
#include "tools/random_functions.h"
#include "tools/timer.h"
#include "io/graph_io.h"
#include "spac/spac.h"
#include "tools/quality_metrics.h"

static void execute_kahip(graph_access &G, PartitionConfig &config);

int main(int argn, char **argv) {
    PartitionConfig partition_config;
    SpacConfig spac_config;
    std::string graph_filename;

    if (parse_spac_parameters(argn, argv, partition_config, spac_config, graph_filename)) {
        return 0;
    }

    std::cout << "graph: " << graph_filename << "\n"
              << "infinity edge weight: " << spac_config.infinity << "\n"
              << "seed: " << partition_config.seed << "\n"
              << "k: " << partition_config.k << std::endl;

    timer t;

    // load input graph
    t.restart();
    graph_access input_graph;
    if (graph_io::readGraphWeighted(input_graph, graph_filename)) {
        return 1;
    }
    std::cout << "input IO took " << t.elapsed() << "\n"
              << "n(input): " << input_graph.number_of_nodes() << "\n"
              << "m(input): " << input_graph.number_of_edges() << std::endl;

    // construct split graph
    t.restart();
    spac splitter(input_graph, spac_config.infinity);
    graph_access &split_graph = splitter.construct_split_graph();
    std::cout << "split graph construction took " << t.elapsed() << "\n"
              << "n(split): " << split_graph.number_of_nodes() << "\n"
              << "m(split): " << split_graph.number_of_edges() << std::endl;

    // partition split graph
    t.restart();
    execute_kahip(split_graph, partition_config);
    std::cout << "kahip took " << t.elapsed() << "\n"
              << "edge cut: " << quality_metrics().edge_cut(split_graph) << std::endl;

    // evaluate edge partition
    t.restart();
    splitter.fix_cut_dominant_edges();
    std::vector <PartitionID> edge_partition = splitter.project_partition();
    unsigned vertex_cut = splitter.calculate_vertex_cut(edge_partition);
    std::cout << "evaluation took " << t.elapsed() << "\n"
              << "vertex cut: " << vertex_cut << std::endl;
    return 0;
}

static void execute_kahip(graph_access &G, PartitionConfig &config) {
    G.set_partition_count(config.k);

    balance_configuration bc;
    bc.configurate_balance(config, G);

    std::srand(static_cast<unsigned int>(config.seed));
    random_functions::setSeed(config.seed);

    graph_partitioner().perform_partitioning(config, G);
}