/******************************************************************************
 * graph_communication.cpp 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
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

#include <mpi.h>

#include "graph_communication.h"

graph_communication::graph_communication() {
                
}

graph_communication::~graph_communication() {
                
}

void graph_communication::broadcast_graph( graph_access & G, unsigned root) {
        int rank       = MPI::COMM_WORLD.Get_rank();

        //first B-Cast number of nodes and number of edges 
        unsigned number_of_nodes = 0;
        unsigned number_of_edges = 0;
 
        std::vector< int > buffer(2,0);
        if(rank == (int)root) {
               buffer[0] = G.number_of_nodes();
               buffer[1] = G.number_of_edges();
        }
        MPI::COMM_WORLD.Bcast(&buffer[0], 2, MPI_INT, root);

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

        MPI::COMM_WORLD.Bcast(xadj,   number_of_nodes+1, MPI_INT, root);
        MPI::COMM_WORLD.Bcast(adjncy, number_of_edges  , MPI_INT, root);
        MPI::COMM_WORLD.Bcast(vwgt,   number_of_nodes  , MPI_INT, root);
        MPI::COMM_WORLD.Bcast(adjwgt, number_of_edges  , MPI_INT, root);

        G.build_from_metis_weighted( number_of_nodes, xadj, adjncy, vwgt, adjwgt); 

        delete[] xadj;
        delete[] adjncy;
        delete[] vwgt;
        delete[] adjwgt;
 
}
