/******************************************************************************
 * parallel_graph_io.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
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

                //static int readGraphWeightedMETIS(parallel_graph_access & G, 
                                //std::string filename, 
                                //PEID peID, 
                                //PEID comm_size, MPI_Comm communicator = MPI_COMM_WORLD);

                static int readGraphWeightedFlexible(parallel_graph_access & G, std::string filename, PEID peID, PEID comm_size, MPI_Comm communicator = MPI_COMM_WORLD); 

                //static int readGraphWeightedMETISFast(parallel_graph_access & G, 
                                //std::string filename, 
                                //PEID peID, 
                                //PEID comm_size, MPI_Comm communicator = MPI_COMM_WORLD);

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

                static int readGraphWeightedMETIS_fixed(parallel_graph_access & G, std::string filename, PEID peID, PEID comm_size, MPI_Comm communicator = MPI_COMM_WORLD); 


};


#endif /* end of include guard: PARALLEL_GRAPH_IO_8HHCKD13 */
