/******************************************************************************
 * label_propagation.cpp 
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
#include <sstream>
#include <math.h>
#include <regex.h>
#include <stdio.h>
#include <string.h> 

#include "data_structure/graph_access.h"
#include "graph_io.h"
#include "parse_parameters.h"
#include "partition/partition_config.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"
#include "coarsening/clustering/size_constraint_label_propagation.h"

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

        timer t;
        graph_io::readGraphWeighted(G, graph_filename);
        std::cout << "io time: " << t.elapsed()  << std::endl;
       
        std::cout <<  "graph has " <<  G.number_of_nodes() <<  " nodes and " <<  G.number_of_edges() <<  " edges"  << std::endl;
        if( partition_config.cluster_upperbound == std::numeric_limits< NodeWeight >::max()/2 ) {
                std::cout <<  "no size-constrained specified" << std::endl;
        } else {
                std::cout <<  "size-constrained set to " <<  partition_config.cluster_upperbound << std::endl;
        }
        
        partition_config.upper_bound_partition = partition_config.cluster_upperbound+1;
        partition_config.cluster_coarsening_factor = 1;
        partition_config.k = 1;
        srand(partition_config.seed);
        random_functions::setSeed(partition_config.seed);

        t.restart();
        // ***************************** perform clustering ***************************************       
        NodeID no_blocks = 0;
        std::vector< NodeID > cluster_id(G.number_of_nodes());
        size_constraint_label_propagation sclp;
        sclp.label_propagation( partition_config, G, cluster_id, no_blocks);
        // ******************************* done clustering  *****************************************       
        std::cout <<  "time spent for partitioning " << t.elapsed()  << std::endl;

        // output some information about the partition that we have computed 
        quality_metrics qm;
        forall_nodes(G, node) {
                G.setPartitionIndex(node, cluster_id[node]);
        } endfor
        
        G.set_partition_count(no_blocks);
        std::cout << "number of clusters/blocks  " << no_blocks << std::endl;
        std::cout << "number of edges between clusters " << qm.edge_cut(G)                 << std::endl;

        // write the clustering  to the disc 
        std::stringstream filename;
        if(!partition_config.filename_output.compare("")) {
                // no output filename given
                filename << "tmpclustering";
        } else {
                filename << partition_config.filename_output;
        }
        graph_io::writeVector(cluster_id, filename.str());

}
