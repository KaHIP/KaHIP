/******************************************************************************
 * buffoon.cpp 
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

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string.h> 
#include <argtable2.h>
#include <regex.h>
#include <mpi.h>
#include <math.h>

#include "buffoon/buffoon.h"
#include "data_structure/graph_access.h"
#include "graph_io.h"
#include "macros_assertions.h"
#include "mpi_tools.h"
#include "parse_parameters.h"
#include "partition/graph_partitioner.h"
#include "partition/partition_config.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"

int main(int argn, char **argv) {


        MPI::Init(argn, argv);    /* starts MPI */

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
        int rank = MPI::COMM_WORLD.Get_rank();
        if(rank == ROOT) {
                graph_io::readGraphWeighted(G, graph_filename);
        }
       
        G.set_partition_count(partition_config.k); 
        NodeWeight largest_graph_weight = 0;
        if( rank == ROOT) {
                forall_nodes(G, node) {
                        largest_graph_weight += G.getNodeWeight(node);
                } endfor
        }
        MPI::COMM_WORLD.Bcast(&largest_graph_weight,         1, MPI_INT, ROOT);
        
        double epsilon = (partition_config.imbalance)/100.0;
        partition_config.largest_graph_weight       = largest_graph_weight;
        partition_config.upper_bound_partition      = (1+epsilon)*ceil(largest_graph_weight/(double)partition_config.k);
        std::cout <<  "upper bound " <<  partition_config.upper_bound_partition  << std::endl;
        partition_config.graph_allready_partitioned = false;
        partition_config.kway_adaptive_limits_beta  = log(largest_graph_weight);

        partition_config.buffoon                    = true;

        if(partition_config.input_partition != "") {
                std::cout <<  "reading input partition" << std::endl;
                graph_io::readPartition(G, partition_config.input_partition);
                partition_config.graph_allready_partitioned = true;
        }

        mpi_tools::non_active_wait_for_root();
        MPI::COMM_WORLD.Bcast(&partition_config.upper_bound_partition,        1, MPI_INT, ROOT);
        MPI::COMM_WORLD.Bcast(&partition_config.largest_graph_weight,         1, MPI_INT, ROOT);
        MPI::COMM_WORLD.Bcast(&partition_config.graph_allready_partitioned,   1, MPI_INT, ROOT);
        MPI::COMM_WORLD.Bcast(&partition_config.kway_adaptive_limits_beta,    1, MPI_INT, ROOT);

        srand(partition_config.seed);
        random_functions::setSeed(partition_config.seed);

        if( rank == ROOT ) {
                std::cout << "log> graph has " <<  G.number_of_nodes() <<  
                             " nodes and " <<  G.number_of_edges() <<  " edges"  << std::endl;
                std::cout << "log> io time: " << t.elapsed()  << std::endl;
        }

        // ***************************** perform partitioning ***************************************       
        t.restart();
        buffoon partitioner;
        partitioner.perform_partitioning(partition_config, G);
        
        // ******************************* done partitioning *****************************************       
        ofs.close();
        std::cout.rdbuf(backup);
        if( rank == ROOT ) {
                std::cout <<  "time spent for partitioning " << t.elapsed()  << std::endl;

                // output some information about the partition that we have computed 
                quality_metrics qm;
                std::cout << "cut \t\t"         << qm.edge_cut(G)       << std::endl;
                std::cout << "finalobjective  " << qm.edge_cut(G)       << std::endl;
                std::cout << "bnd \t\t"         << qm.boundary_nodes(G) << std::endl;
                std::cout << "balance \t"       << qm.balance(G)        << std::endl;
                std::cout << "finalbalance \t"  << qm.balance(G)        << std::endl;


                // write the partition to the disc 
                std::string partition("tmppartition");
                std::stringstream noparts;
                noparts << "tmppartition" << partition_config.k;
                graph_io::writePartition(G, noparts.str());
        }

        MPI::Finalize();
        
}
