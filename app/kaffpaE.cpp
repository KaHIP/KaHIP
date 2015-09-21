/******************************************************************************
 * kaffpaE.cpp 
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
#include <mpi.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h> 

#include "algorithms/cycle_search.h"
#include "balance_configuration.h"
#include "data_structure/graph_access.h"
#include "graph_io.h"
#include "macros_assertions.h"
#include "parallel_mh/parallel_mh_async.h"
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

        partition_config.LogDump(stdout);
        partition_config.graph_filename = graph_filename.substr( graph_filename.find_last_of( '/' ) +1 );

        graph_access G;     

        timer t;
        graph_io::readGraphWeighted(G, graph_filename);

        std::cout << "io time: " << t.elapsed()  << std::endl;

        omp_set_num_threads(1);
        G.set_partition_count(partition_config.k); 
        partition_config.kaffpaE = true; // necessary for balance configuration 

        balance_configuration bc;
        bc.configurate_balance( partition_config, G);

        std::vector<PartitionID> input_partition;
        if(partition_config.input_partition != "") {
                std::cout <<  "reading input partition" << std::endl;
                graph_io::readPartition(G, partition_config.input_partition);
                partition_config.graph_allready_partitioned = true;

                input_partition.resize(G.number_of_nodes());

                forall_nodes(G, node) {
                        input_partition[node] = G.getPartitionIndex(node);
                } endfor
        }

        t.restart();

        parallel_mh_async mh;
        mh.perform_partitioning(partition_config, G);

        
        int rank = MPI::COMM_WORLD.Get_rank();
        if( rank == ROOT ) {
                std::cout <<  "time spent for partitioning " << t.elapsed()  << std::endl;
                std::cout <<  "time spent in neg. cycle detection " <<  cycle_search::total_time  << std::endl;
                std::cout <<  "time spent in neg. cycle detection (rel) " <<  (cycle_search::total_time/t.elapsed()*100)  << std::endl;

                // output some information about the partition that we have computed 
                quality_metrics qm;
                EdgeWeight cut =  qm.edge_cut(G);
                std::cout << "cut \t\t"         << cut                            << std::endl;
                std::cout << "finalobjective  " << cut                            << std::endl;
                std::cout << "bnd \t\t"         << qm.boundary_nodes(G)           << std::endl;
                std::cout << "balance \t"       << qm.balance(G)                  << std::endl;
                std::cout << "max_comm_vol \t"  << qm.max_communication_volume(G) << std::endl;

                // write the partition to the disc 
                std::stringstream filename;
                if(!partition_config.filename_output.compare("")) {
                        // no output filename given
                        filename << "tmppartition" << partition_config.k;
                } else {
                        filename << partition_config.filename_output;
                }

                graph_io::writePartition(G, filename.str());
        }
 
        MPI::Finalize();
}
