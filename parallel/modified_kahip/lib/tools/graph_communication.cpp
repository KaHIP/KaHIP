/******************************************************************************
 * graph_communication.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <mpi.h>

#include "graph_communication.h"

graph_communication::graph_communication() {
                
}

graph_communication::~graph_communication() {
                
}

void graph_communication::broadcast_graph( graph_access & G, unsigned root) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        //first B-Cast number of nodes and number of edges 
        unsigned number_of_nodes = 0;
        unsigned number_of_edges = 0;
 
        std::vector< int > buffer(2,0);
        if(rank == (int)root) {
               buffer[0] = G.number_of_nodes();
               buffer[1] = G.number_of_edges();
        }
        MPI_Bcast(&buffer[0], 2, MPI_INT, root, MPI_COMM_WORLD);

        number_of_nodes = buffer[0];
        number_of_edges = buffer[1];

        int* xadj;        
        int* adjncy;
        int* vwgt;        
        int* adjwgt;

        if( rank == (int)root) {
                xadj           = G.UNSAFE_metis_style_xadj_array();
                adjncy         = G.UNSAFE_metis_style_adjncy_array();

                vwgt           = G.UNSAFE_metis_style_vwgt_array();
                adjwgt         = G.UNSAFE_metis_style_adjwgt_array();
        } else {
                xadj   = new int[number_of_nodes+1];
                adjncy = new int[number_of_edges];

                vwgt   = new int[number_of_nodes];
                adjwgt = new int[number_of_edges];
        }

        MPI_Bcast(xadj,   number_of_nodes+1, MPI_INT, root, MPI_COMM_WORLD);
        MPI_Bcast(adjncy, number_of_edges  , MPI_INT, root, MPI_COMM_WORLD);
        MPI_Bcast(vwgt,   number_of_nodes  , MPI_INT, root, MPI_COMM_WORLD);
        MPI_Bcast(adjwgt, number_of_edges  , MPI_INT, root, MPI_COMM_WORLD);

        G.build_from_metis_weighted( number_of_nodes, xadj, adjncy, vwgt, adjwgt); 

        delete[] xadj;
        delete[] adjncy;
        delete[] vwgt;
        delete[] adjwgt;
 
}
