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
#include "partition/uncoarsening/separator/area_bfs.h"
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
       
        partition_config.k = 2;
        G.set_partition_count(partition_config.k); 
 
        std::cout <<  "imbalance is set to " <<  partition_config.imbalance << "%" << std::endl;
        balance_configuration bc;
        bc.configurate_balance( partition_config, G);

        srand(partition_config.seed);
        random_functions::setSeed(partition_config.seed);

        std::cout <<  "graph has " <<  G.number_of_nodes() <<  " nodes and " <<  G.number_of_edges() <<  " edges"  << std::endl;
        // ***************************** perform partitioning ***************************************       
        t.restart();
        graph_partitioner partitioner;
        quality_metrics qm;
        area_bfs::m_deepth.resize(G.number_of_nodes());
        forall_nodes(G, node) {
                area_bfs::m_deepth[node] = 0;
        } endfor
        

        std::cout <<  "computing a node separator"  << std::endl;
        partition_config.mode_node_separators = true;
        partitioner.perform_partitioning(partition_config, G);

        // ******************************* done partitioning *****************************************       
        ofs.close();
        std::cout.rdbuf(backup);
        std::cout <<  "time spent to compute vertex separator " << t.elapsed()  << std::endl;
       
        NodeWeight ns_size = 0;
        forall_nodes(G, node) {
                if(G.getPartitionIndex(node) == G.getSeparatorBlock()) {
                        ns_size++;
                }
        } endfor
        G.set_partition_count(3);
        std::cout << "balance "        << qm.balance_separator(G) << std::endl;
        std::cout << "separator size " << ns_size << std::endl;

        std::stringstream filename;
        if(!partition_config.filename_output.compare("")) {
                filename << "tmpseparator" << partition_config.k;
        } else {
                filename << partition_config.filename_output;
        }

        graph_io::writePartition(G, filename.str());

#ifdef NDEBUG
        // check wether it is still a separator
        bool not_a_separator = false;
        forall_nodes(G, node) {
                if( G.getPartitionIndex(node) == 0 ) {
                        forall_out_edges(G, e, node) {
                                NodeID target = G.getEdgeTarget(e);
                                if( G.getPartitionIndex(target) == 1 ) {
					not_a_separator = true;
					break;
                                }
                        } endfor
                }
        } endfor

	if( not_a_separator ) {
		std::cout <<  "not a separator -- please report this bug"  << std::endl;
        }
#endif
}
