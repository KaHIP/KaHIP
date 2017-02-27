/******************************************************************************
 * parallel_graph_io.h
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

#ifndef PARALLEL_GRAPH_IO_8HHCKD13
#define PARALLEL_GRAPH_IO_8HHCKD13

#include <fstream>
#include <iostream>

#include "definitions.h"
#include "data_structure/parallel_graph_access.h"
#include "partition_config.h"

#define MAXLINE	50000000
class parallel_graph_io {
        public:
                parallel_graph_io();
                virtual ~parallel_graph_io();

                static int readGraphWeighted(PPartitionConfig & config, parallel_graph_access & G, 
                                std::string filename, 
                                PEID peID, 
                                PEID comm_size, MPI_Comm communicator = MPI_COMM_WORLD);

                static int readGraphWeightedMETIS(parallel_graph_access & G, 
                                std::string filename, 
                                PEID peID, 
                                PEID comm_size, MPI_Comm communicator = MPI_COMM_WORLD);

                static int readGraphWeightedFlexible(parallel_graph_access & G, std::string filename, PEID peID, PEID comm_size, MPI_Comm communicator = MPI_COMM_WORLD); 

                static int readGraphWeightedMETISFast(parallel_graph_access & G, 
                                std::string filename, 
                                PEID peID, 
                                PEID comm_size, MPI_Comm communicator = MPI_COMM_WORLD);

                static int readGraphBinary(PPartitionConfig & config, parallel_graph_access & G, 
                                std::string filename, 
                                PEID peID, 
                                PEID comm_size, MPI_Comm communicator = MPI_COMM_WORLD);

                static int writeGraphParallelSimple(parallel_graph_access & G, 
                                std::string filename, MPI_Comm communicator = MPI_COMM_WORLD);

                static int writeGraphWeightedParallelSimple(parallel_graph_access & G, 
                                std::string filename, MPI_Comm communicator = MPI_COMM_WORLD);

                static int writeGraphWeightedSequentially(complete_graph_access & G, 
                                std::string filename);

                static int writeGraphSequentially(complete_graph_access & G, 
                                std::string filename);

                static int writeGraphSequentially(complete_graph_access & G, 
                                std::ofstream & f);

                static int writeGraphSequentiallyBinary(complete_graph_access & G, std::string filename);

                static int writeGraphExternallyBinary(std::string intput_filename, std::string output_filename);

};


#endif /* end of include guard: PARALLEL_GRAPH_IO_8HHCKD13 */
