/******************************************************************************
 * evaluator.cpp 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
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
#include <regex.h>
#include <string.h> 

#include "data_structure/graph_access.h"
#include "graph_io.h"
#include "macros_assertions.h"
#include "parse_parameters.h"
#include "partition/partition_config.h"
#include "quality_metrics.h"

int main(int argn, char **argv) {

        PartitionConfig partition_config;
        std::string graph_filename;

        bool is_graph_weighted = false;
        bool suppress_output   = false;
        bool recursive         = false;
       
        int ret_code = parse_parameters(argn, argv, 
                                        partition_config, 
                                        graph_filename, 
                                        is_graph_weighted, 
                                        suppress_output, recursive); 

        if(ret_code) {
                return 0;
        }

        graph_access G;     
        graph_io::readGraphWeighted(G, graph_filename);

        G.set_partition_count(partition_config.k); 
 
        std::vector<PartitionID> input_partition;
        if(partition_config.input_partition != "") {
                std::cout <<  "reading input partition" << std::endl;
                graph_io::readPartition(G, partition_config.input_partition);
        } else {
                std::cout <<  "Please specify an input partition using the --input_partition flag."  << std::endl;
                exit(0);
        }

        std::cout <<  "graph has " <<  G.number_of_nodes() <<  " nodes and " <<  G.number_of_edges() <<  " edges"  << std::endl;
        quality_metrics qm;

        std::cout << "cut \t\t"         << qm.edge_cut(G)                 << std::endl;
        std::cout << "no boundary vertices \t\t" << qm.boundary_nodes(G)           << std::endl;
        std::cout << "balance \t"       << qm.balance(G)                  << std::endl;
        std::cout << "balance based on edges \t"       << qm.balance_edges(G)                  << std::endl;
        std::cout << "max comm vol \t"  << qm.max_communication_volume(G) << std::endl;
        std::cout << "min comm vol \t"  << qm.min_communication_volume(G) << std::endl;
        std::cout << "total comm vol \t"  << qm.total_communication_volume(G) << std::endl;
}
