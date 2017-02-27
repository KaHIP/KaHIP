/******************************************************************************
 * toolbox.cpp
 *
 * Source of KaHIP -- Karlsruhe High Quality Graph Partitioning 
 ******************************************************************************
 * Copyright (C) 2017 Christian Schulz 
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
#include <iomanip>
#include <mpi.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h> 

#include "communication/mpi_tools.h"
#include "communication/dummy_operations.h"
#include "data_structure/parallel_graph_access.h"
#include "distributed_partitioning/distributed_partitioner.h"
#include "io/parallel_graph_io.h"
#include "io/parallel_vector_io.h"
#include "macros_assertions.h"
#include "parse_parameters.h"
#include "partition_config.h"
#include "random_functions.h"
#include "timer.h"
#include "tools/distributed_quality_metrics.h"

int main(int argn, char **argv) {

        MPI_Init(&argn, &argv);    /* starts MPI */

        PPartitionConfig partition_config;
        std::string graph_filename;

        int ret_code = parse_parameters(argn, argv, 
                        partition_config, 
                        graph_filename); 

        if(ret_code) {
                MPI_Finalize();
                return 0;
        }

        int rank, size;
        MPI_Comm communicator = MPI_COMM_WORLD; 
        MPI_Comm_rank( communicator, &rank);
        MPI_Comm_size( communicator, &size);

        MPI_Comm_rank( communicator, &rank);
        MPI_Comm_size( communicator, &size);

        partition_config.stop_factor /= partition_config.k;
        if(rank != 0) partition_config.seed = partition_config.seed*size+rank; 

        srand(partition_config.seed);

        parallel_graph_access G(communicator);
        parallel_graph_io::readGraphWeighted(partition_config, G, graph_filename, rank, size, communicator);
        parallel_vector_io pvio;
        pvio.readPartition(partition_config, G, partition_config.input_partition_filename);

        G.printMemoryUsage(std::cout);

        MPI_Barrier(communicator);

        if(partition_config.converter_evaluate) {
                distributed_quality_metrics qm;
                EdgeWeight edge_cut = qm.edge_cut( G, communicator );
                double balance  = qm.balance( partition_config, G, communicator );
                double balance_load  = qm.balance_load( partition_config, G, communicator );
                double balance_load_dist  = qm.balance_load_dist( partition_config, G, communicator );

                if( rank == ROOT ) {
                        std::cout << "log>" << "=====================================" << std::endl;
                        std::cout << "log>" << "============Evaluation Result========" << std::endl;
                        std::cout << "log>" << "=====================================" << std::endl;
                        std::cout <<  "log>final edge cut " <<  edge_cut  << std::endl;
                        std::cout <<  "log>final balance "  <<  balance   << std::endl;
                        std::cout <<  "log>final balance load "  <<  balance_load   << std::endl;
                        std::cout <<  "log>final balance load dist "  <<  balance_load_dist   << std::endl;
                }
                qm.comm_vol( partition_config, G, communicator );
        }

        if( partition_config.save_partition ) {
                if(rank == ROOT) std::cout <<  "saving text partition"  << std::endl;
                parallel_vector_io pvio;
                std::string filename("tmppartition.txtp");
                pvio.writePartitionSimpleParallel(G, filename);
        }

        if( partition_config.save_partition_binary ) {
                if(rank == ROOT) std::cout <<  "saving binary partition"  << std::endl;
                parallel_vector_io pvio;
                std::string filename("tmppartition.binp");
                pvio.writePartitionBinaryParallelPosix(partition_config, G, filename);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
}
