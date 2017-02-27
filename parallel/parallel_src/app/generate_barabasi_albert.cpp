/******************************************************************************
 * generate_barabasi_albert.cpp
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
#include <mpi.h>

#include "communication/mpi_tools.h"
#include "communication/dummy_operations.h"
#include "data_structure/parallel_graph_access.h"
#include "io/parallel_graph_io.h"
#include "io/parallel_vector_io.h"
#include "macros_assertions.h"
#include "parse_parameters.h"
#include "partition_config.h"
#include "timer.h"
#include "graph_generation/generate_barabasi_albert.h"

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
        MPI_Comm_rank( MPI_COMM_WORLD, &rank);
        MPI_Comm_size( MPI_COMM_WORLD, &size);

        timer t;
        MPI_Barrier(MPI_COMM_WORLD);
        {
                t.restart();
                if( rank == ROOT ) std::cout <<  "running collective dummy operations ";
                dummy_operations dop;
                dop.run_collective_dummy_operations();
        }
        MPI_Barrier(MPI_COMM_WORLD);

        if( rank == ROOT ) {
                std::cout <<  "took " <<  t.elapsed()  << std::endl;
        }
        t.restart();


        for( int i = 0; i < 10; i++) {
                parallel_graph_access G;
                generate_barabasialbert gba;

                if( partition_config.compute_degree_sequence_k_first ) {
                        gba.generate_k_deghist(partition_config, G, partition_config.n > 0, graph_filename.compare(""));
                        break;
                } else {
                        if(partition_config.generate_ba_32bit) {
                                gba.generate_32bit(partition_config, G, partition_config.n > 0, graph_filename.compare(""));
                        } else {
                                gba.generate(partition_config, G, partition_config.n > 0, graph_filename.compare(""));
                        }
                        if( graph_filename.compare("") && i == 9) {
                                if( rank == ROOT ) { std::cout <<  "writing to disk"  << std::endl; }
                                parallel_graph_io::writeGraphParallelSimple( G, graph_filename );
                        }
                }

                MPI_Barrier(MPI_COMM_WORLD);

        }

        
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();

}
