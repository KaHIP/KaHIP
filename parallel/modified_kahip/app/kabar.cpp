/******************************************************************************
 * kabar.cpp 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013 Christian Schulz <christian.schulz@kit.edu>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
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

#include <algorithm>
#include <argtable2.h>
#include <iostream>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h> 

#include "algorithms/cycle_search.h"
#include "data_structure/graph_access.h"
#include "graph_io.h"
#include "macros_assertions.h"
#include "parse_parameters.h"
#include "partition/graph_partitioner.h"
#include "partition/partition_config.h"
#include "partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.h"
#include "partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_commons.h"
#include "partition_snapshooter.h"
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
                                        partition_config, graph_filename, 
                                        is_graph_weighted, suppress_output, 
                                        recursive); 

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

        NodeWeight largest_graph_weight = 0;
        forall_nodes(G, node) {
                largest_graph_weight += G.getNodeWeight(node);
        } endfor

        double epsilon                              = partition_config.imbalance/100.0;
        partition_config.upper_bound_partition      = (1+epsilon)*ceil(largest_graph_weight/(double)partition_config.k);
        partition_config.largest_graph_weight       = largest_graph_weight;
        partition_config.graph_allready_partitioned = false;
        partition_config.kway_adaptive_limits_beta  = log(largest_graph_weight);


        if(partition_config.input_partition != "") {
                std::cout <<  "reading input partition" << std::endl;
                graph_io::readPartition(G, partition_config.input_partition);
                partition_config.graph_allready_partitioned = true;
        } else {
                std::cout <<  "please specify an input partition."  << std::endl;
                exit(0);
        }

        srand(partition_config.seed);
        random_functions::setSeed(partition_config.seed);

        std::cout <<  "graph has " <<  G.number_of_nodes() <<  " nodes and " <<  G.number_of_edges() <<  " edges"  << std::endl;
        // ***************************** perform tests ***************************************       
        quality_metrics qm;

        
        t.restart();
        complete_boundary boundary(&G);
        boundary.build();
        cycle_refinement cr;
        //now obtain the quotient graph
        Gain overall_gain = cr.perform_refinement(partition_config, G, boundary);

        std::cout <<  "overall gain " <<  overall_gain  << std::endl;
        // ******************************** done *********************************************       
        ofs.close();
        std::cout.rdbuf(backup);

        std::cout << "time spent "                               << t.elapsed()                                << std::endl;
        std::cout << "time spent in neg. cycle detection "       << cycle_search::total_time                   << std::endl;
        std::cout << "time spent in neg. cycle detection (rel) " << (cycle_search::total_time/t.elapsed()*100) << std::endl;

        // output some information about the partition that we have computed 
        std::cout << "cut \t\t"         << qm.edge_cut(G)       << std::endl;
        std::cout << "finalobjective  " << qm.edge_cut(G)       << std::endl;
        std::cout << "bnd \t\t"         << qm.boundary_nodes(G) << std::endl;
        std::cout << "balance \t"       << qm.balance(G)        << std::endl;
        std::cout << "finalbalance \t"  << qm.balance(G)        << std::endl;

        std::cout <<  "conflicts: "<<  advanced_models::conflicts  << std::endl;
        partition_snapshooter* pss = partition_snapshooter::getInstance();
        pss->flush_buffer();

        // write the partition to the disc 
        std::stringstream noparts;
        noparts << "tmppartition_bal" << partition_config.k;
        graph_io::writePartition(G, noparts.str());


}
