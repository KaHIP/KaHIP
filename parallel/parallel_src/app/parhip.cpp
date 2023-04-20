/******************************************************************************
 * parhip.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <argtable3.h>
#include <iostream>
#include <math.h>
#include <iomanip>
#include <mpi.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h> 

#include "communication/mpi_tools.h"
#include "communication/dummy_operations.h"
#include "data_structure/parallel_graph_access.h"
#include "distributed_partitioning/distributed_partitioner.h"
#include "io/parallel_graph_io.h"
#include "io/parallel_vector_io.h"
#include "macros_assertions.h"
#include "parse_parameters.h"
#include "partition_config.h"
#include "random_functions.h"
#include "timer.h"
#include "tools/distributed_quality_metrics.h"

int main(int argn, char **argv) {

        MPI_Init(&argn, &argv);    /* starts MPI */

        PPartitionConfig partition_config;
        std::string graph_filename;

        int ret_code = parse_parameters(argn, argv, 
                        partition_config, 
                        graph_filename); 

        if(ret_code) {
                MPI_Finalize();
                return 0;
        }

        int rank, size;
        MPI_Comm communicator = MPI_COMM_WORLD; 
        MPI_Comm_rank( communicator, &rank);
        MPI_Comm_size( communicator, &size);

        timer t;
        MPI_Barrier(MPI_COMM_WORLD);
        {
                t.restart();
                if( rank == ROOT ) std::cout <<  "running collective dummy operations ";
                dummy_operations dop;
                dop.run_collective_dummy_operations();
        }
        MPI_Barrier(MPI_COMM_WORLD);

        if( rank == ROOT ) {
                std::cout <<  "took " <<  t.elapsed()  << std::endl;
        }

        if( communicator != MPI_COMM_NULL) {
                MPI_Comm_rank( communicator, &rank);
                MPI_Comm_size( communicator, &size);

                if(rank == ROOT) {
                        PRINT(std::cout <<  "log> cluster coarsening factor is set to " <<  partition_config.cluster_coarsening_factor  << std::endl;)
                }

                partition_config.stop_factor /= partition_config.k;
                if(rank != 0) partition_config.seed = partition_config.seed*size+rank; 

                srand(partition_config.seed);

                parallel_graph_access G(communicator);
                parallel_graph_io::readGraphWeighted(partition_config, G, graph_filename, rank, size, communicator);
                //parallel_graph_io::readGraphWeightedFlexible(G, graph_filename, rank, size, communicator);
                if( rank == ROOT ) std::cout <<  "took " <<  t.elapsed()  << std::endl;
                if( rank == ROOT ) std::cout <<  "n:" <<  G.number_of_global_nodes() << " m: " <<  G.number_of_global_edges()  << std::endl;

                random_functions::setSeed(partition_config.seed);
                parallel_graph_access::set_comm_rounds( partition_config.comm_rounds/size );
                parallel_graph_access::set_comm_rounds_up( partition_config.comm_rounds/size);
                distributed_partitioner::generate_random_choices( partition_config );

                G.printMemoryUsage(std::cout);

                //compute some stats
                EdgeWeight interPEedges = 0;
                EdgeWeight localEdges = 0;
                NodeWeight localWeight = 0;
                forall_local_nodes(G, node) {
                        localWeight += G.getNodeWeight(node);
                        forall_out_edges(G, e, node) {
                                NodeID target = G.getEdgeTarget(e);
                                if(!G.is_local_node(target)) {
                                        interPEedges++;
                                } else {
                                        localEdges++;
                                }
                        } endfor
                } endfor

                EdgeWeight globalInterEdges = 0;
                EdgeWeight globalIntraEdges = 0;
                EdgeWeight globalWeight = 0;
                MPI_Reduce(&interPEedges, &globalInterEdges, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOT, communicator);
                MPI_Reduce(&localEdges, &globalIntraEdges, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOT, communicator);
                MPI_Allreduce(&localWeight, &globalWeight, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);

                if( rank == ROOT ) {
                        std::cout <<  "log> ghost edges " <<  globalInterEdges/(double)G.number_of_global_edges() << std::endl;
                        std::cout <<  "log> local edges " <<  globalIntraEdges/(double)G.number_of_global_edges() << std::endl;
                }

                t.restart();
                double epsilon = (partition_config.inbalance)/100.0;
                if( partition_config.vertex_degree_weights ) {
                        NodeWeight total_load = G.number_of_global_edges()+G.number_of_global_edges();
                        partition_config.number_of_overall_nodes = G.number_of_global_nodes();
                        partition_config.upper_bound_partition   = (1+epsilon)*ceil(total_load/(double)partition_config.k);

                        forall_local_nodes(G, node) {
                                G.setNodeWeight(node, G.getNodeDegree(node)+1);
                        } endfor

                } else {
                        partition_config.number_of_overall_nodes = G.number_of_global_nodes();
                        partition_config.upper_bound_partition   = (1+epsilon)*ceil(globalWeight/(double)partition_config.k);
                        if( rank == ROOT) {
                                std::cout <<  "upper bound on blocks " << partition_config.upper_bound_partition  << std::endl;
                        }
                }


                distributed_partitioner dpart;
                dpart.perform_partitioning( communicator, partition_config, G);

                MPI_Barrier(communicator);

                double running_time = t.elapsed();
                distributed_quality_metrics qm;
                EdgeWeight edge_cut = qm.edge_cut( G, communicator );
                double balance  = qm.balance( partition_config, G, communicator );
                PRINT(double balance_load  = qm.balance_load( partition_config, G, communicator );)
                PRINT(double balance_load_dist  = qm.balance_load_dist( partition_config, G, communicator );)

                if( rank == ROOT ) {
                        std::cout << "log>" << "=====================================" << std::endl;
                        std::cout << "log>" << "============AND WE R DONE============" << std::endl;
                        std::cout << "log>" << "=====================================" << std::endl;
                        std::cout <<  "log>total partitioning time elapsed " <<  running_time << std::endl;
                        std::cout <<  "log>final edge cut " <<  edge_cut  << std::endl;
                        std::cout <<  "log>final balance "  <<  balance   << std::endl;
                        PRINT(std::cout <<  "log>final balance load "  <<  balance_load   << std::endl;)
                        PRINT(std::cout <<  "log>final balance load dist "  <<  balance_load_dist   << std::endl;)
                }
                PRINT(qm.comm_vol( partition_config, G, communicator );)
                PRINT(qm.comm_vol_dist( G, communicator );)


#ifndef NDEBUG
                MPI_Status st; int flag; 
                MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, communicator, &flag, &st);
                while( flag ) {
                        std::cout <<  "attention: still incoming messages! rank " <<  rank <<  " from " <<  st.MPI_SOURCE << std::endl;
                        int message_length;
                        MPI_Get_count(&st, MPI_UNSIGNED_LONG_LONG, &message_length);
                        MPI_Status rst;
                        std::vector<NodeID> message; message.resize(message_length);
                        MPI_Recv( &message[0], message_length, MPI_UNSIGNED_LONG_LONG, st.MPI_SOURCE, st.MPI_TAG, communicator, &rst); 
                        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, communicator, &flag, &st);
                };
#endif

                if( partition_config.save_partition ) {
                        parallel_vector_io pvio;
                        std::string filename("tmppartition.txtp");
                        pvio.writePartitionSimpleParallel(G, filename);
                }

                if( partition_config.save_partition_binary ) {
                        parallel_vector_io pvio;
                        std::string filename("tmppartition.binp");
                        pvio.writePartitionBinaryParallelPosix(partition_config, G, filename);
                }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
}
