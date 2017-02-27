/******************************************************************************
 * generate_rgg.h
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


#ifndef GENERATE_RGG_PLHS3WMW
#define GENERATE_RGG_PLHS3WMW

#include "data_structure/parallel_graph_access.h"
#include "error_codes.h"
#include "generate_grid.h"
#include "rgg/generate_rgg.h"
#include "graph_io.h"
#include "macros_assertions.h"
#include "static_distributed_graph.h"
#include "static_local_graph.h"
#include "static_local_graph_properties.h"
#include "util_rand.h"

using namespace distributed_graph;

// This struct stores the program's options.
struct ProgramOptions
{

  // Number of vertices and edges for the R-MAT algorithm.
  long n;

  // The radius for the unit disk graph generator.
  double radius;

  // The seed to use for pseudo random number generation.
  int seed;
};

class generator_rgg {
public:
        generator_rgg() ;
        virtual ~generator_rgg() ;

        void generate(PPartitionConfig & config, parallel_graph_access & G) {
                ProgramOptions options;
                options.n = pow(2, config.log_num_verts);
                options.radius = 0.55*sqrt(log(options.n)/(double)options.n);
                options.seed = config.seed;

                std::vector< std::vector< NodeID > > messages;
                NodeID nmbNodes = options.n;
                NodeID nmbEdges;
                int global_rank = MPI::COMM_WORLD.Get_rank();

                // create a communicator with a square number of processors 
                MPI::Intracomm  comm_world, comm_worker;
                MPI::Group group_world, group_worker;

                comm_world = MPI::COMM_WORLD;
                group_world = comm_world.Get_group();

                int size = MPI::COMM_WORLD.Get_size();
                int sq = floor(sqrt(size)); sq *= sq;

                std::vector< int > excluded_ranks;
                for( int i = sq; i < size; i++) {
                        excluded_ranks.push_back(i);
                }

                int num_excluded = excluded_ranks.size();
                group_worker = group_world.Excl(num_excluded, &excluded_ranks[0]);  //[> process 0 not member <]
                comm_worker = comm_world.Create(group_worker);
                
                if(global_rank < sq) {
                        StaticDistributedGraph *graph;
                        RandomGeometricGraphGenerator generator;
                        RandomGeometricGraphGeneratorConfig config(RGG_UNIT_SQUARE, options.n, options.radius);

                        generator.Init(&comm_worker, config);
                        std::tr1::mt19937 mt;
                        SeedMersenneTwisterDistributedly(&comm_worker, options.seed, &mt);

                        generator.Run(&mt, &graph);
                        if( global_rank == 0 ) { 
                                std::cout << "n: " <<  graph->global_n()  << " m: " <<  graph->global_m()  << std::endl;
                        }
                        nmbEdges = graph->global_m();
                        const VertexId *vertex_splitters = graph->vertex_splitters();

                        // Get a shortcut to the current static local graph.
                        const StaticLocalGraph *local_graph = graph->local_graph();
                        UnsafeStaticLocalGraphAccess local_access(local_graph);

                        // now build up messages for all PEs
                        messages.resize(size);

                        NodeID divisor  = ceil( graph->global_n() / (double)size);
                        //std::cout <<  "divisor " <<  divisor  << std::endl;
                        int rank = graph->communicator()->Get_rank();
                        for( long i = 0, iend = local_graph->n(); i < iend; i++) {
                                VertexId u = i+vertex_splitters[rank];
                                PEID peID = u / divisor;
                                
                                messages[peID].push_back(u);
                                messages[peID].push_back(local_graph->out_degree(u));
                                for( long j = 0, jend = local_graph->out_degree(u); j < jend; j++) {
                                        VertexId v = local_graph->neighbour(u,j);
                                        messages[peID].push_back(v);
                                }
                        }

                        delete graph;
                        for( PEID peID = 0; peID < (PEID)size; peID++) {
                                if( messages[peID].size() == 0 ){
                                        messages[peID].push_back(std::numeric_limits<NodeID>::max());
                                }

                                MPI::COMM_WORLD.Isend( &messages[peID][0], messages[peID].size(), MPI_UNSIGNED_LONG_LONG, peID, 0);
                        }

                }

                MPI::COMM_WORLD.Bcast(&nmbEdges,  1, MPI_UNSIGNED_LONG_LONG, 0);

                PEID peID = global_rank;
                ULONG from = peID     * ceil(nmbNodes / (double)size);
                ULONG to   = (peID+1) * ceil(nmbNodes / (double)size) - 1;
                to         = std::min(to, nmbNodes-1);
                ULONG local_no_nodes = to - from + 1;

                std::vector< std::vector< NodeID > > local_edge_lists;
                local_edge_lists.resize(local_no_nodes);

                std::vector< std::vector< NodeID > >  inc_messages;
                inc_messages.resize(size);

                PEID counter = 0;
                EdgeID edge_counter = 0;

                while( counter < (PEID)(sq )) {
                        // wait for incomming message of an adjacent processor
                        MPI::Status st;
                        while(MPI::COMM_WORLD.Iprobe(MPI::ANY_SOURCE,0,st)) {
                                ULONG message_length = st.Get_count(MPI_UNSIGNED_LONG_LONG);
                                std::vector<NodeID> incmessage; incmessage.resize(message_length);
                                //std::cout <<  "receiving " <<  message_length  << std::endl;

                                MPI::COMM_WORLD.Recv( &incmessage[0], message_length, MPI_UNSIGNED_LONG_LONG, st.Get_source(), 0); 
                                counter++;

                                // now integrate the changes
                                if( incmessage[0] == std::numeric_limits< NodeID >::max()) continue; // nothing to do

                                for( ULONG i = 0; i < message_length; ) {
                                        NodeID source = incmessage[i];
                                        NodeID degree = incmessage[i+1];
                                        if( degree > 0 ) {
                                                NodeID start = i+2;
                                                NodeID end = i+2 + degree;
                                                if( source < from ) {
                                                        std::cout <<  "node does not belong here"  << std::endl;
                                                        exit(0);
                                                }
                                                for( ULONG j = start; j < end; j++) {
                                                        local_edge_lists[source-from].push_back(incmessage[j]);
                                                        edge_counter++;
                                                }
                                                i = end;
                                        } else {
                                                i+=2;
                                        }
                                        
                                }
                        }
                }

                MPI::COMM_WORLD.Barrier();
                messages.resize(0);
                std::vector< std::vector< NodeID > >(messages).swap(messages);

                std::cout <<  "done generating and communicating -- now building"  << std::endl;
                G.start_construction(local_no_nodes, 2*edge_counter, nmbNodes, nmbEdges);
                G.set_range(from, to);

                std::vector< NodeID > vertex_dist( size+1, 0 );
                for( PEID peID = 0; peID <= size; peID++) {
                        vertex_dist[peID] = std::min(nmbNodes, (NodeID) (peID * ceil(nmbNodes / (double)size))); // from positions
                }
                G.set_range_array(vertex_dist);

                for (NodeID i = 0; i < local_no_nodes; ++i) {
                        NodeID node = G.new_node();
                        G.setNodeWeight(node, 1);
                        G.setNodeLabel(node, from+node);
                        G.setSecondPartitionIndex(node, 0);

                        for( ULONG j = 0; j < local_edge_lists[i].size(); j++) {
                                NodeID target = local_edge_lists[i][j]; 
                                EdgeID e = G.new_edge(node, target);
                                G.setEdgeWeight(e, 1);
                        }
                }

                G.finish_construction();
        }
};


#endif /* end of include guard: GENERATE_RGG_PLHS3WMW */
