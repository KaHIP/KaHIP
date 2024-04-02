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
        NodeID number_of_nodes = 0;
        EdgeID number_of_edges = 0;

        MPI_Datatype MNodeID = (sizeof(NodeID) == 4) ? MPI_INT : MPI_LONG_LONG;
        MPI_Datatype MEdgeID = (sizeof(EdgeID) == 4) ? MPI_INT : MPI_LONG_LONG;

        if(rank == (int)root) {
               number_of_nodes = G.number_of_nodes();
               number_of_edges = G.number_of_edges();
        }

        MPI_Bcast(&number_of_nodes, 1, MNodeID, root, MPI_COMM_WORLD);
        MPI_Bcast(&number_of_edges, 1, MEdgeID, root, MPI_COMM_WORLD);

        EdgeID* xadj;
        int* adjncy;
        int* vwgt;        
        int* adjwgt;

        if( rank == (int)root) {
                xadj           = G.UNSAFE_metis_style_xadj_array();
                adjncy         = G.UNSAFE_metis_style_adjncy_array();

                vwgt           = G.UNSAFE_metis_style_vwgt_array();
                adjwgt         = G.UNSAFE_metis_style_adjwgt_array();
        } else {
                xadj   = new EdgeID[number_of_nodes+1];
                adjncy = new int[number_of_edges];

                vwgt   = new int[number_of_nodes];
                adjwgt = new int[number_of_edges];
        }

        // FIXME: When 64-bit need to use Bcast_c
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
