/******************************************************************************
 * generate_kronecker.cpp 
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

#include <vector>
#include <algorithm>
#include "generate_kronecker.h"
#include "generate_g500_edgelist.h"
#include "helpers.h"

generate_kronecker::generate_kronecker() {
                
}

generate_kronecker::~generate_kronecker() {
                
}

void generate_kronecker::generate_kronecker_graph( PPartitionConfig & config, parallel_graph_access & G) {
        std::vector< source_target_pair > edge_list;
        generate_g500_edgelist generator;
        generator.generate_kronecker_edgelist( config, edge_list );

        if(config.kronecker_internal_only) return;
        ULONG my_edges     = edge_list.size();
        ULONG global_edges = 0;
        MPI_Allreduce(&my_edges, &global_edges, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        PEID rank, size;
        MPI_Comm_rank( MPI_COMM_WORLD, &rank);
        MPI_Comm_size( MPI_COMM_WORLD, &size);
        
        if(rank == ROOT) {
                std::cout <<  "remove self loops. insert backward edges: now " <<  global_edges  <<  " edges " << std::endl;
        }

        // now send each processor its part of the graph
        // first build the messages
        NodeID global_n = INT64_C(1) << config.log_num_verts;
        NodeID divisor  = ceil( global_n / (double)size);

        helpers helper;
        helper.filter_duplicates(edge_list, 
                        [](const source_target_pair & lhs, const source_target_pair& rhs) -> bool { 
                        return (lhs.source <  rhs.source || (lhs.source == rhs.source &&  lhs.target <  rhs.target)); 
                        }, 
                        [](const source_target_pair & lhs, const source_target_pair& rhs) -> bool { 
                        return (lhs.source ==  rhs.source &&  lhs.target ==  rhs.target); 
                        }); // can be omitted if memory is not an issue


        std::vector< std::vector< NodeID > >  messages;
        messages.resize(size);
        for( ULONG i = 0; i < edge_list.size(); i++) {
                PEID peID = edge_list[i].source / divisor;
                messages[ peID ].push_back( edge_list[i].source );
                messages[ peID ].push_back( edge_list[i].target );
        }

        edge_list.resize(0);
        std::vector< source_target_pair >(edge_list).swap(edge_list);

        // now flood the network
        // TODO: this should become an all to all operation here 
        for( PEID peID = 0; peID < size; peID++) {
                if( peID != rank ) {
                        if( messages[peID].size() == 0 ){
                                // length 1 encode no message
                                messages[peID].push_back(0);
                        }

                        MPI_Request rq; 

                        MPI_Isend( &messages[peID][0], 
                                        messages[peID].size(), 
                                        MPI_UNSIGNED_LONG_LONG, 
                                        peID, peID, MPI_COMM_WORLD, &rq);
                }
        }

        std::vector< source_target_pair > final_local_edge_list;
        if( messages[rank].size() != 0 ) {
                for( ULONG i = 0; i < messages[rank].size()-1; i+=2) {
                        source_target_pair st;
                        st.source = messages[rank][i];
                        st.target = messages[rank][i+1];
                        final_local_edge_list.push_back( st );
                }
                messages[rank].resize(0);
                std::vector< NodeID >(messages[rank]).swap(messages[rank]);
        }

        PEID counter = 0;
        while( counter < size - 1) {
                // wait for incomming message of an adjacent processor
                MPI_Status st;
                MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
                
                int message_length;
                MPI_Get_count(&st, MPI_UNSIGNED_LONG_LONG, &message_length);
                std::vector<NodeID> incmessage; incmessage.resize(message_length);

                MPI_Status rst;
                MPI_Recv( &incmessage[0], message_length, MPI_UNSIGNED_LONG_LONG, st.MPI_SOURCE, rank, MPI_COMM_WORLD, &rst); 
                
                counter++;

                // now integrate the changes
                if(message_length == 1) continue; // nothing to do

                for( int i = 0; i < message_length-1; i+=2) {
                        source_target_pair st;
                        st.source = incmessage[i];
                        st.target = incmessage[i+1];
                        final_local_edge_list.push_back(st);
                }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        
        //improve memory footprint
        for( PEID peID = 0; peID < size; peID++) {
                messages[peID].resize(0);
                std::vector< NodeID >(messages[peID]).swap(messages[peID]);
        }

        helper.filter_duplicates(final_local_edge_list, 
                        [](const source_target_pair & lhs, const source_target_pair& rhs) -> bool { 
                           return (lhs.source <  rhs.source || (lhs.source == rhs.source &&  lhs.target <=  rhs.target)); 
                        }, 
                        [](const source_target_pair & lhs, const source_target_pair& rhs) -> bool { 
                           return (lhs.source ==  rhs.source &&  lhs.target ==  rhs.target); 
                        });

        // create the graph -- meanwhile filter parallel edges
        ULONG from           = rank     * ceil(global_n / (double)size);
        ULONG to             = (rank+1) * ceil(global_n / (double)size) - 1;
        to = std::min(to, global_n -1);

        ULONG local_no_nodes = to - from + 1;
        G.start_construction(local_no_nodes, final_local_edge_list.size(), global_n, 2*global_edges);
        G.set_range(from, to);

        std::vector< NodeID > vertex_dist( size+1, 0 );
        for( PEID peID = 0; peID <= size; peID++) {
                vertex_dist[peID] = std::min(global_n, (NodeID) (peID * ceil(global_n/ (double)size))); // from positions
        }
        G.set_range_array(vertex_dist);

        EdgeID pos = 0;
        for (NodeID i = 0; i < local_no_nodes; ++i) {
                NodeID node = G.new_node();
                NodeID globalID = from+node;
                G.setNodeWeight(node, 1);
                G.setNodeLabel(node, globalID);
               
                if( final_local_edge_list.size() == 0 ) {continue;}
                while( final_local_edge_list[pos].source == globalID) {
                        // we still have to create edges
                        NodeID target = final_local_edge_list[pos].target;
                        EdgeID e = G.new_edge(node, target);
                        G.setEdgeWeight(e, 1);
               
                        //last_target = target; 
                        pos++;
                        if( pos == final_local_edge_list.size()) break;
                }
        }

        G.finish_construction();

        // reset number global edges since we have removed parallel edges
        my_edges     = G.number_of_local_edges();
        global_edges = 0;

        MPI_Allreduce(&my_edges, &global_edges, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD); 
        
        G.set_number_of_global_edges( global_edges );

        // and we are done
        if( rank == ROOT ) std::cout <<  "graph with " <<  global_edges  << " global edges has been successfully created and distributed"  << std::endl;

        MPI_Barrier(MPI_COMM_WORLD);
}
