/******************************************************************************
 * mpi_tools.cpp
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

#include <mpi.h>
#include <unistd.h>

#include "io/parallel_vector_io.h"
#include "io/parallel_graph_io.h"
#include "mpi_tools.h"

mpi_tools::mpi_tools() {
                
}

mpi_tools::~mpi_tools() {
                

}

// currently this method is for debugging purposses only
// later on this may be a parallel io routine
void mpi_tools::collect_and_write_labels( MPI_Comm communicator, PPartitionConfig & config, parallel_graph_access & G) {
        int rank, size;
        MPI_Comm_rank( communicator, &rank);
        MPI_Comm_size( communicator, &size);

        std::vector< NodeID > labels;

        if( rank == ROOT ) {
                labels.resize(G.number_of_global_nodes());
                forall_local_nodes(G, node) {
                        labels[node] = G.getNodeLabel(node);
                } endfor
        } else {
                //pack the data
                forall_local_nodes(G, node) {
                        labels.push_back(G.getGlobalID(node));
                        labels.push_back(G.getNodeLabel(node));
                } endfor
        }
        
        if( rank == ROOT ) {
                int counter = 0;
                while( counter < size-1) {
                        // wait for incomming message of an adjacent processor
                        int flag; MPI_Status st;
                        MPI_Iprobe(MPI_ANY_SOURCE, rank, communicator, &flag, &st);
                         
                        while( flag ) {
                                int message_length;
                                MPI_Get_count(&st, MPI_UNSIGNED_LONG_LONG, &message_length);
                                std::vector<NodeID> message; message.resize(message_length);

                                MPI_Status rst;
                                MPI_Recv( &message[0], message_length, MPI_UNSIGNED_LONG_LONG, st.MPI_SOURCE, st.MPI_TAG, communicator, &rst); 
                                counter++;

                                for( int i = 0; i < message_length-1; i+=2) {
                                        NodeID global_id = message[i];
                                        NodeID label     = message[i+1];

                                        labels[global_id] = label;
                                }
                                MPI_Iprobe(MPI_ANY_SOURCE, rank, communicator, &flag, &st);
                        }
                }

        } else {
                MPI_Request rq;
                MPI_Isend( &labels[0], labels.size(), MPI_UNSIGNED_LONG_LONG, ROOT, rank+12*size, communicator, &rq);
        }

        if( rank == ROOT ) {
                std::string clustering_filename("tmpclustering");
                parallel_vector_io pvio;
                pvio.writeVectorSequentially(labels, clustering_filename);
        }
        MPI_Barrier(communicator);
}


void mpi_tools::collect_parallel_graph_to_local_graph( MPI_Comm communicator, PPartitionConfig & config, 
                                                       parallel_graph_access & G,
                                                       complete_graph_access & Q) {

        int rank, size;
        MPI_Comm_rank( communicator, &rank);
        MPI_Comm_size( communicator, &size);
        
        std::vector< NodeID > message;

        if( rank == ROOT ) {
                Q.start_construction( G.number_of_global_nodes(), G.number_of_global_edges(), 
                                      G.number_of_global_nodes(), G.number_of_global_edges(), false); // no update of comm_rounds!
                Q.set_range(0, G.number_of_global_nodes()); // this graph should contain all global edges
                forall_local_nodes(G, node) {
                        NodeID cur_node = Q.new_node();
                        Q.setNodeWeight(cur_node, G.getNodeWeight(node));
                        Q.setSecondPartitionIndex(cur_node, G.getSecondPartitionIndex(node));
                        forall_out_edges(G, e, node) {
                                NodeID target = G.getEdgeTarget(e);
                                EdgeID e_bar = Q.new_edge(cur_node, G.getGlobalID(target));
                                Q.setEdgeWeight(e_bar, G.getEdgeWeight(e));
                        } endfor
                } endfor
        } else {
                // layout: no local  nodes, no local edges, node_1, pidx, weight, degree,its edges: e_1, w_1, ...,node_2, ... 
                message.push_back(G.number_of_local_nodes());
                forall_local_nodes(G, node) {
                        //message.push_back(G.getGlobalID(node));
                        message.push_back(G.getSecondPartitionIndex(node));
                        message.push_back(G.getNodeWeight(node));
                        message.push_back(G.getNodeDegree(node));
                        forall_out_edges(G, e, node) {
                                NodeID target = G.getEdgeTarget(e);
                                message.push_back(G.getGlobalID(target));
                                message.push_back(G.getEdgeWeight(e));
                        } endfor
                } endfor
        }
        
        if( rank == ROOT) {
                for( int i = 1; i < size; i++) {
                        int flag; MPI_Status st;
                        MPI_Iprobe(i, 13*size, communicator, &flag, &st);
                        
                        while(!flag) { MPI_Iprobe(i, 13*size, communicator, &flag, &st); }
                                
                        int message_length;
                        MPI_Get_count(&st, MPI_UNSIGNED_LONG_LONG, &message_length);
                        std::vector<NodeID> rmessage; rmessage.resize(message_length);

                        MPI_Status rst;
                        MPI_Recv( &rmessage[0], message_length, MPI_UNSIGNED_LONG_LONG, st.MPI_SOURCE, st.MPI_TAG, communicator, &rst); 
                        
                        NodeID no_nodes = rmessage[0];
                        NodeID pos = 1;
                        for( ULONG node = 0; node < no_nodes; node++) {
                                NodeID cur_node = Q.new_node();
                                Q.setSecondPartitionIndex(cur_node, rmessage[pos++]);
                                Q.setNodeWeight(cur_node, rmessage[pos++]);

                                EdgeID degree = rmessage[pos++];
                                for( ULONG e = 0; e < degree; e++) {
                                        EdgeID e_bar = Q.new_edge(cur_node, rmessage[pos++]);
                                        Q.setEdgeWeight(e_bar, rmessage[pos++]);
                                }
                        }

                }
        } else {
                MPI_Request rq; 
                MPI_Isend( &message[0], message.size(), MPI_UNSIGNED_LONG_LONG, ROOT, 13*size, communicator, &rq);
        }

        if( rank == ROOT ) {
                Q.finish_construction();
        }

        MPI_Barrier(communicator);
}



void mpi_tools::distribute_local_graph( MPI_Comm communicator, PPartitionConfig & config, 
                                        complete_graph_access & G) {

        int rank;
        MPI_Comm_rank( communicator, &rank);

        //first B-Cast number of nodes and number of edges 
        ULONG number_of_nodes = 0;
        ULONG number_of_edges = 0;

        std::vector< int > buffer(2,0);
        if(rank == (int)ROOT) {
                buffer[0] = G.number_of_global_nodes();
                buffer[1] = G.number_of_global_edges();
        }
        MPI_Bcast(&buffer[0], 2, MPI_INT, ROOT, communicator);

        number_of_nodes = buffer[0];
        number_of_edges = buffer[1];

        int* xadj;        
        int* adjncy;
        int* vwgt;        
        int* adjwgt;

        if( rank == (int)ROOT) {
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
        MPI_Bcast(xadj, number_of_nodes+1, MPI_INT, ROOT, communicator);
        MPI_Bcast(adjncy, number_of_edges, MPI_INT, ROOT, communicator);
        MPI_Bcast(vwgt, number_of_nodes, MPI_INT, ROOT, communicator);
        MPI_Bcast(adjwgt, number_of_edges, MPI_INT, ROOT, communicator);

        G.build_from_metis_weighted( number_of_nodes, xadj, adjncy, vwgt, adjwgt); 

        delete[] xadj;
        delete[] adjncy;
        delete[] vwgt;
        delete[] adjwgt;
}

void mpi_tools::alltoallv( void * sendbuf, 
                ULONG sendcounts[], ULONG displs[], 
                const MPI_Datatype & sendtype, void * recvbuf,
                ULONG recvcounts[], ULONG rdispls[],
                const MPI_Datatype & recvtype, MPI_Comm communicator ) {

        int rank, size;
        MPI_Comm_rank( communicator, &rank);
        MPI_Comm_size( communicator, &size);
        
        bool no_special_case = true;
        for( int i = 0; i < size && no_special_case; i++) {
                if( sendcounts[i] > std::numeric_limits< int >::max()) no_special_case = false;
                if( recvcounts[i] > std::numeric_limits< int >::max()) no_special_case = false;
        }
        if( displs[size]  > std::numeric_limits< int >::max()) no_special_case = false;
        if( rdispls[size] > std::numeric_limits< int >::max()) no_special_case = false;

        if( no_special_case ) {
                int sbktsize[size];
                int rbktsize[size];
                int sdispl[size+1];
                int rdispl[size+1];

                for( int i = 0; i < size; i++) {
                        sbktsize[i] = sendcounts[i];
                        rbktsize[i] = recvcounts[i];
                }

                for( int i = 0; i <= size; i++) {
                        sdispl[i] = displs[i];
                        rdispl[i] = rdispls[i];
                }

                MPI_Alltoallv(sendbuf, sbktsize, sdispl, MPI_UNSIGNED_LONG_LONG, 
                              recvbuf, rbktsize, rdispl, MPI_UNSIGNED_LONG_LONG, communicator);
        } else {
                if( rank == ROOT ) { std::cout <<  "special case all to all with counts > sizeof(int)! not tested yet!"  << std::endl; exit(0);}
        }
}

