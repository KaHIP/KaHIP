/******************************************************************************
 * node_separator.cpp 
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
#include <iostream>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h> 

#include "balance_configuration.h"
#include "data_structure/graph_access.h"
#include "graph_io.h"
#include "macros_assertions.h"
#include "parse_parameters.h"
#include "partition/graph_partitioner.h"
#include "partition/partition_config.h"
#include "partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.h"
#include "partition/uncoarsening/separator/vertex_separator_algorithm.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"

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

        std::streambuf* backup = std::cout.rdbuf();
        std::ofstream ofs;
        ofs.open("/dev/null");
        if(suppress_output) {
                std::cout.rdbuf(ofs.rdbuf()); 
        }

        partition_config.LogDump(stdout);
        graph_access G;     

        timer t;
        graph_io::readGraphWeighted(G, graph_filename);
        std::cout << "io time: " << t.elapsed()  << std::endl;
       
        G.set_partition_count(partition_config.k); 
 
        balance_configuration bc;
        bc.configurate_balance( partition_config, G);

        srand(partition_config.seed);
        random_functions::setSeed(partition_config.seed);

        std::cout <<  "graph has " <<  G.number_of_nodes() <<  " nodes and " <<  G.number_of_edges() <<  " edges"  << std::endl;
        // ***************************** perform partitioning ***************************************       
        t.restart();
        graph_partitioner partitioner;
        quality_metrics qm;

        std::cout <<  "computing a node separator"  << std::endl;

        partition_config.mode_node_separators = true;
        partition_config.use_fullmultigrid    = false;
        partition_config.use_wcycles          = false;
        partition_config.matching_type        = MATCHING_GPA;
        partitioner.perform_partitioning(partition_config, G);

        // ******************************* done partitioning *****************************************       
        ofs.close();
        std::cout.rdbuf(backup);
        std::cout <<  "time spent to compute vertex separator " << t.elapsed()  << std::endl;
       
        // output some information about the partition that we have computed 
        //std::cout << "initial separator size " << separator.size()        << std::endl;
        //std::cout << "new separator size "     << output_separator.size() << std::endl;
        //std::cout << "cut \t\t"                << qm.edge_cut(G)          << std::endl;
        //std::cout << "finalobjective  "        << qm.edge_cut(G)          << std::endl;

        
}
