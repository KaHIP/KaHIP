/******************************************************************************
 * readbgf.cpp
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

#include <stdio.h>
#include <iostream>
#include "io/parallel_graph_io.h"
#include "partition_config.h"
#include "configuration.h"

using namespace std;

const long fileTypeVersionNumber = 3;
const long header_count    = 3;

int main(int argn, char **argv)
{

        MPI_Init(&argn, &argv);    /* starts MPI */

        int rank, size;
        MPI_Comm_rank( MPI_COMM_WORLD, &rank);
        MPI_Comm_size( MPI_COMM_WORLD, &size);

        if(argn != 2) {
                if( rank == ROOT ) {
                        std::cout <<  "usage: " ;
                        std::cout <<  "readbgf bfg_file"  << std::endl;
                }
                MPI_Finalize();
                return 0;
        }


        if( rank == ROOT ) {
                std::cout <<  "program reads a BGF (binary graph format) file and prints it into dummy. "  << std::endl;
        }
        string filename(argv[1]);

        configuration cfg;
        PPartitionConfig config;
        cfg.standard(config);

        parallel_graph_access G;
        parallel_graph_io pgio;
        pgio.readGraphBinary(config, G, filename, rank, size);

        string output_filename("dummy");
        parallel_graph_io::writeGraphParallelSimple(G, output_filename);

        MPI_Finalize();
        return 0;
}

