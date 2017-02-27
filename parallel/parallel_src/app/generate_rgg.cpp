/******************************************************************************
 * parallel_label_compress.cpp 
 *
 * Source of the Parallel Partitioning Program
 ******************************************************************************
 * Copyright (C) 2014 Christian Schulz <christian.schulz@kit.edu>
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

#include <argtable2.h>
#include <mpi.h>

#include "communication/mpi_tools.h"
#include "data_structure/parallel_graph_access.h"
#include "io/parallel_graph_io.h"
#include "io/parallel_vector_io.h"
#include "macros_assertions.h"
#include "parse_parameters.h"
#include "partition_config.h"
#include "timer.h"
#include "graph_generation/generate_rgg.h"

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
        t.restart();

        parallel_graph_access G;
        generator_rgg grgg;
        grgg.generate(partition_config, G);

        MPI_Barrier(MPI_COMM_WORLD);

        if( rank == ROOT ) {
                std::cout <<  "graph generation took " <<  t.elapsed()  << std::endl;
                std::cout <<  "writing to disk"  << std::endl;
        }
        
        parallel_graph_io::writeGraphParallelSimple( G, graph_filename );

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();

}
