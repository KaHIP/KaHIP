/******************************************************************************
 * kaHIP_interface.cpp
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

#include <iostream>
#include <mpi.h>
#include "kaHIP_interface.h"
#include "../lib/data_structure/graph_access.h"
#include "../lib/io/graph_io.h"
#include "../lib/tools/timer.h"
#include "../lib/tools/quality_metrics.h"
#include "../lib/tools/macros_assertions.h"
#include "../lib/tools/random_functions.h"
#include "../lib/parallel_mh/parallel_mh_async.h"
#include "../lib/partition/partition_config.h"
#include "../lib/partition/graph_partitioner.h"
#include "../lib/partition/uncoarsening/separator/vertex_separator_algorithm.h"
#include "../app/configuration.h"

using namespace std;

void internal_build_graph( PartitionConfig & partition_config, 
                           int* n, 
                           int* vwgt, 
                           int* xadj, 
                           int* adjcwgt, 
                           int* adjncy,
                           graph_access & G) {
        G.build_from_metis(*n, xadj, adjncy); 
        G.set_partition_count(partition_config.k); 
 
        srand(partition_config.seed);
        random_functions::setSeed(partition_config.seed);
       
        if(vwgt != NULL) {
                forall_nodes(G, node) {
                        G.setNodeWeight(node, vwgt[node]);
                } endfor
        }

        if(adjcwgt != NULL) {
                forall_edges(G, e) {
                        G.setEdgeWeight(e, adjcwgt[e]);
                } endfor 
        }

        partition_config.largest_graph_weight = 0;
        forall_nodes(G, node) {
                partition_config.largest_graph_weight += G.getNodeWeight(node);
        } endfor
        
        double epsilon = partition_config.imbalance/100;
        partition_config.upper_bound_partition = ceil((1+epsilon)*partition_config.largest_graph_weight/(double)partition_config.k);
        partition_config.graph_allready_partitioned = false;

}


void internal_kaffpa_call(PartitionConfig & partition_config, 
                          bool suppress_output, 
                          int* n, 
                          int* vwgt, 
                          int* xadj, 
                          int* adjcwgt, 
                          int* adjncy, 
                          int* nparts, 
                          double* imbalance, 
                          int* edgecut, 
                          int* part) {

        streambuf* backup = cout.rdbuf();
        ofstream ofs;
        ofs.open("/dev/null");
        if(suppress_output) {
               cout.rdbuf(ofs.rdbuf()); 
        }

        partition_config.imbalance = 100*(*imbalance);
        graph_access G;     
        internal_build_graph( partition_config, n, vwgt, xadj, adjcwgt, adjncy, G);

        
        timer t;
        graph_partitioner partitioner;
        partitioner.perform_partitioning(partition_config, G);
        std::cout <<  "partioning took " <<  t.elapsed()  << std::endl;

        forall_nodes(G, node) {
                part[node] = G.getPartitionIndex(node);
        } endfor

        quality_metrics qm;
        *edgecut = qm.edge_cut(G);

        ofs.close();
        cout.rdbuf(backup);
}



void kaffpa(int* n, 
                   int* vwgt, 
                   int* xadj, 
                   int* adjcwgt, 
                   int* adjncy, 
                   int* nparts, 
                   double* imbalance, 
                   bool suppress_output, 
                   int seed,
                   int mode,
                   int* edgecut, 
                   int* part) {
        configuration cfg;
        PartitionConfig partition_config;
        partition_config.k = *nparts;

        switch( mode ) {
                case FAST: 
                        cfg.fast(partition_config);
                        break;
                case ECO: 
                        cfg.eco(partition_config);
                        break;
                case STRONG: 
                        cfg.strong(partition_config);
                        break;
                case FASTSOCIAL: 
                        cfg.fastsocial(partition_config);
                        break;
                case ECOSOCIAL: 
                        cfg.ecosocial(partition_config);
                        break;
                case STRONGSOCIAL: 
                        cfg.strongsocial(partition_config);
                        break;
                default: 

                        cfg.eco(partition_config);
                        break;
        }

        partition_config.seed = seed;
        internal_kaffpa_call(partition_config, suppress_output, n, vwgt, xadj, adjcwgt, adjncy, nparts, imbalance, edgecut, part);
}


void internal_nodeseparator_call(PartitionConfig & partition_config, 
                          bool suppress_output, 
                          int* n, 
                          int* vwgt, 
                          int* xadj, 
                          int* adjcwgt, 
                          int* adjncy, 
                          int* nparts, 
                          double* imbalance, 
                          int* num_nodeseparator_vertices, 
                          int** separator) {

        //first perform std partitioning using KaFFPa
        streambuf* backup = cout.rdbuf();
        ofstream ofs;
        ofs.open("/dev/null");
        if(suppress_output) {
               cout.rdbuf(ofs.rdbuf()); 
        }

        partition_config.imbalance = 100*(*imbalance);
        graph_access G;     
        internal_build_graph( partition_config, n, vwgt, xadj, adjcwgt, adjncy, G);

        
        graph_partitioner partitioner;
        partitioner.perform_partitioning(partition_config, G);

        // now compute a node separator from the partition of the graph
        complete_boundary boundary(&G);
        boundary.build();

        vertex_separator_algorithm vsa;
        std::vector<NodeID> internal_separator;
        vsa.compute_vertex_separator(partition_config, G, boundary, internal_separator);

        // copy to output variables
        *num_nodeseparator_vertices =  internal_separator.size();
        *separator = new int[*num_nodeseparator_vertices];
        for( unsigned int i = 0; i < internal_separator.size(); i++) {
                (*separator)[i] = internal_separator[i];
        }

        ofs.close();
        cout.rdbuf(backup);
}


void node_separator(int* n, 
                          int* vwgt, 
                          int* xadj, 
                          int* adjcwgt, 
                          int* adjncy, 
                          int* nparts, 
                          double* imbalance, 
                          bool suppress_output, 
                          int seed,
                          int mode,
                          int* num_separator_vertices, 
                          int** separator) {
        configuration cfg;
        PartitionConfig partition_config;
        partition_config.k = *nparts;

        switch( mode ) {
                case FAST: 
                        cfg.fast(partition_config);
                        break;
                case ECO: 
                        cfg.eco(partition_config);
                        break;
                case STRONG: 
                        cfg.strong(partition_config);
                        break;
                case FASTSOCIAL: 
                        cfg.fastsocial(partition_config);
                        break;
                case ECOSOCIAL: 
                        cfg.ecosocial(partition_config);
                        break;
                case STRONGSOCIAL: 
                        cfg.strongsocial(partition_config);
                        break;
                default: 
                        cfg.eco(partition_config);
                        break;
        }
        partition_config.seed = seed;

        internal_nodeseparator_call(partition_config, suppress_output, n, vwgt, xadj, adjcwgt, adjncy, nparts, imbalance, num_separator_vertices, separator);
}


void kaffpaE(int* n, 
                   int* vwgt, 
                   int* xadj, 
                   int* adjcwgt, 
                   int* adjncy, 
                   int* nparts, 
                   double* imbalance, 
                   bool suppress_output, 
                   bool graph_partitioned, 
                   int time_limit,
                   int seed,
                   int mode, // 0 == strong, 1 == eco, 2 == fast
                   MPI_Comm communicator, 
                   int* edgecut, 
                   double* balance,
                   int* part) {

        configuration cfg;
        PartitionConfig partition_config;
        partition_config.k = *nparts;
        cfg.standard(partition_config);

        switch( mode ) {
                case FAST: 
                        cfg.fast(partition_config);
                        break;
                case ECO: 
                        cfg.eco(partition_config);
                        break;
                case STRONG: 
                        cfg.strong(partition_config);
                        break;
                case FASTSOCIAL: 
                        cfg.fastsocial(partition_config);
                        break;
                case ULTRAFASTSOCIAL: 
                        cfg.fastsocial(partition_config);
			partition_config.ultra_fast_kaffpaE_interfacecall = true;
                        break;
                case ECOSOCIAL: 
                        cfg.ecosocial(partition_config);
                        break;
                case STRONGSOCIAL: 
                        cfg.strongsocial(partition_config);
                        break;
                default: 
                        cfg.eco(partition_config);
                        break;
        }

        partition_config.seed      = seed;
        partition_config.k         = *nparts;
        partition_config.imbalance = 100*(*imbalance);
        partition_config.time_limit = time_limit;
        partition_config.kabapE = false;

        graph_access G;     
        internal_build_graph( partition_config, n, vwgt, xadj, adjcwgt, adjncy, G);

        partition_config.kway_adaptive_limits_beta  = log(partition_config.largest_graph_weight);

        if(graph_partitioned) {
                forall_nodes(G, node) {
                         G.setPartitionIndex(node, part[node]);
                } endfor
        }
        
        partition_config.graph_allready_partitioned  = graph_partitioned;
        partition_config.no_new_initial_partitioning = graph_partitioned;

        parallel_mh_async mh(communicator);
        mh.perform_partitioning(partition_config, G);

        forall_nodes(G, node) {
                part[node] = G.getPartitionIndex(node);
        } endfor

        quality_metrics qm;
        *edgecut = qm.edge_cut(G);
        *balance = qm.balance(G);

        //ofs.close();
        //cout.rdbuf(backup);
}


