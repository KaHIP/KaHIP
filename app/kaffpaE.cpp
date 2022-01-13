/******************************************************************************
 * kaffpaE.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <argtable3.h>
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
#include <omp.h>

int main(int argn, char **argv) {

        MPI_Init(&argn, &argv);    /* starts MPI */

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
        if( partition_config.imbalance < 1 ) {
                partition_config.kabapE = true;
        }

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

        
        int rank, size;
        MPI_Comm communicator = MPI_COMM_WORLD; 
        MPI_Comm_rank( communicator, &rank);
        MPI_Comm_size( communicator, &size);

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
 
        MPI_Finalize();
}
