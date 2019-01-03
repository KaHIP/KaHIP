/******************************************************************************
 * graph2binary_external.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <stdio.h>
#include <iostream>
#include "io/parallel_graph_io.h"

using namespace std;

const long fileTypeVersionNumber = 2;
const long header_count    = 3;

int main(int argn, char **argv)
{
        std::cout <<  "program converts a METIS graph file into a binary (distributed graph format) file. "  << std::endl;

        MPI_Init(&argn, &argv);    /* starts MPI */

        int rank, size;
        MPI_Comm_rank( MPI_COMM_WORLD, &rank);
        MPI_Comm_size( MPI_COMM_WORLD, &size);

        if(argn != 3) {
                if( rank == ROOT ) {
                        std::cout <<  "usage: " ;
                        std::cout <<  "graph2binary_external metisfile outputfilename"  << std::endl;
                }
                MPI_Finalize();
                return 0;
        }

        if( size > 1 ) {
                std::cout <<  "currently only one process supported."  << std::endl;
                MPI_Finalize();
                return 0;
        }

        string graph_filename(argv[1]);
        string filename(argv[2]);

        std::cout <<  "Reading and writing graph " << graph_filename  << std::endl;
        parallel_graph_io::writeGraphExternallyBinary(graph_filename, filename);
        
        MPI_Finalize();
        return 0;
}

