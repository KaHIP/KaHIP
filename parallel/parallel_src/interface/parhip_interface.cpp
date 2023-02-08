#include <iostream>

#include "parhip_interface.h"
#include "parallel_graph_io.h"
#include "configuration.h"
#include "distributed_partitioning/distributed_partitioner.h"
#include "tools/distributed_quality_metrics.h"
#include "random_functions.h"


// 3% imbalance should be specified as imbalance = 0.03 
void ParHIPPartitionKWay(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt,
                         int *nparts, double* imbalance, bool suppress_output, int seed, int mode, int *edgecut, idxtype *part, 
                         MPI_Comm *comm) {



        std::streambuf* backup = std::cout.rdbuf();
        std::ofstream ofs;
        ofs.open("/dev/null");

        if(suppress_output) {
                std::cout.rdbuf(ofs.rdbuf()); 
        }

        PEID rank, size;
        MPI_Comm_rank( *comm, &rank);
        MPI_Comm_size( *comm, &size);

        //building internal graph data structure
        idxtype local_number_of_nodes = vtxdist[rank+1] - vtxdist[rank];
        idxtype local_number_of_edges = xadj[local_number_of_nodes];
        idxtype number_of_nodes = vtxdist[size];

        std::vector< NodeID > vertex_weights(local_number_of_nodes,1);
        NodeWeight local_overall_node_weight = local_number_of_nodes;
        NodeWeight global_node_weight = number_of_nodes;
        if( vwgt != NULL ) {
                local_overall_node_weight = 0;
                global_node_weight = 0;
                for( unsigned long long i = 0; i < local_number_of_nodes; i++) {
                        vertex_weights[i] = vwgt[i];
                        local_overall_node_weight += vwgt[i]; 
                }
                MPI_Allreduce(&local_overall_node_weight, &global_node_weight, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, *comm); 
        }
       
        //// pe p obtains nodes p*ceil(n/size) to (p+1)floor(n/size) and the edges
        idxtype from = vtxdist[rank];
        idxtype to   = vtxdist[rank+1]-1;


        unsigned long long global_number_of_edges = 0;
        MPI_Allreduce(&local_number_of_edges, &global_number_of_edges, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, *comm); 

        parallel_graph_access G(*comm);
        G.start_construction(local_number_of_nodes, local_number_of_edges, number_of_nodes, global_number_of_edges);
        G.set_range(from, to);
        std::vector< NodeID > vertex_dist( size+1, 0 );
        for( PEID peID = 0; peID <= size; peID++) {
                vertex_dist[peID] = vtxdist[peID];
        }
        G.set_range_array(vertex_dist);

        if( adjwgt != NULL ) {
                for (NodeID i = 0; i < local_number_of_nodes; ++i) {
                        NodeID node = G.new_node();
                        G.setNodeWeight(node, vertex_weights[i]);
                        G.setNodeLabel(node, from+node);
                        G.setSecondPartitionIndex(node, 0);

                        for (ULONG j = xadj[i]; j < xadj[i + 1]; j++) {
                                EdgeID e = G.new_edge(node, adjncy[j]);
                                G.setEdgeWeight(e, adjwgt[j]);
                        }
                }
        } else {
                for (NodeID i = 0; i < local_number_of_nodes; ++i) {
                        NodeID node = G.new_node();
                        G.setNodeWeight(node, vertex_weights[i]);
                        G.setNodeLabel(node, from+node);
                        G.setSecondPartitionIndex(node, 0);

                        for (ULONG j = xadj[i]; j < xadj[i + 1]; j++) {
                                EdgeID e = G.new_edge(node, adjncy[j]);
                                G.setEdgeWeight(e, 1);
                        }
                }       
        }

        G.finish_construction();

        PPartitionConfig partition_config;
        configuration cfg;
        cfg.standard(partition_config);

        switch( mode ) {
                case FASTMESH: 
                        cfg.fast(partition_config);
                        partition_config.cluster_coarsening_factor = 20000;
                        break;
                case ULTRAFASTMESH: 
                        cfg.ultrafast(partition_config);
                        partition_config.cluster_coarsening_factor = 20000;
                        break;
                case ECOMESH: 
                        cfg.eco(partition_config);
                        partition_config.cluster_coarsening_factor = 20000;
                        break;
                case FASTSOCIAL: 
                        cfg.fast(partition_config);
                        break;
                case ECOSOCIAL: 
                        cfg.eco(partition_config);
                        break;
                case ULTRAFASTSOCIAL: 
                        cfg.ultrafast(partition_config);
                        break;
                default: 
                        cfg.fast(partition_config);
                        break;
        }

        partition_config.k = *nparts;
        partition_config.seed = seed;
        partition_config.stop_factor /= partition_config.k;
        if(rank != 0) partition_config.seed = partition_config.seed*size+rank; 

        srand(partition_config.seed);

        random_functions::setSeed(partition_config.seed);
        parallel_graph_access::set_comm_rounds( partition_config.comm_rounds/size );
        parallel_graph_access::set_comm_rounds_up( partition_config.comm_rounds/size);
        distributed_partitioner::generate_random_choices( partition_config );

        timer t; 
        partition_config.inbalance = 100*(*imbalance);
        double epsilon = (partition_config.inbalance)/100.0;
        partition_config.number_of_overall_nodes = G.number_of_global_nodes();
        partition_config.upper_bound_partition   = (1+epsilon)*ceil(global_node_weight/(double)partition_config.k);

        distributed_partitioner dpart;
        dpart.perform_partitioning( *comm, partition_config, G);
        MPI_Barrier(*comm);
        double running_time = t.elapsed();

        ofs.close();
        std::cout.rdbuf(backup);
        
        distributed_quality_metrics qm;
        *edgecut = qm.edge_cut( G, *comm );

        if (!suppress_output) {
                double balance  = qm.balance( partition_config, G, *comm );
                if (rank == 0) {
                        std::cout << "log>" << "=====================================" << std::endl;
                        std::cout << "log>" << "============AND WE R DONE============" << std::endl;
                        std::cout << "log>" << "=====================================" << std::endl;
                        std::cout <<  "log>total partitioning time elapsed " << running_time << std::endl;
                        std::cout <<  "log>final edge cut " <<  *edgecut  << std::endl;
                        std::cout <<  "log>final balance "  <<  balance   << std::endl;
                }
        }

        for (NodeID i = 0; i < local_number_of_nodes; ++i) {
                part[i] = G.getNodeLabel(i);
        }
}
