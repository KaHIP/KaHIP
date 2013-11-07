/******************************************************************************
 * kaffpa_interface.cpp 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013 Christian Schulz <christian.schulz@kit.edu>
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
#include "../app/configuration.h"

using namespace std;

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

        graph_access G;     

        G.build_from_metis(*n, xadj, adjncy); 
        G.set_partition_count(*nparts); 
        partition_config.k         = *nparts;
        partition_config.imbalance = 100*(*imbalance);
 
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

        double epsilon                        = partition_config.imbalance;
        partition_config.largest_graph_weight = 0;
        forall_nodes(G, node) {
                partition_config.largest_graph_weight += G.getNodeWeight(node);
        } endfor
        
        partition_config.upper_bound_partition      = ceil((1+epsilon)*partition_config.largest_graph_weight/(double)partition_config.k);
        partition_config.graph_allready_partitioned = false;
        
        cout <<  "performing partitioning"  << endl;
        graph_partitioner partitioner;
        partitioner.perform_partitioning(partition_config, G);

        forall_nodes(G, node) {
                part[node] = G.getPartitionIndex(node);
        } endfor

        quality_metrics qm;
        *edgecut = qm.edge_cut(G);

        ofs.close();
        cout.rdbuf(backup);
}



void kaffpa_strong(int* n, 
                   int* vwgt, 
                   int* xadj, 
                   int* adjcwgt, 
                   int* adjncy, 
                   int* nparts, 
                   double* imbalance, 
                   bool suppress_output, 
                   int seed,
                   int* edgecut, 
                   int* part) {
        configuration cfg;
        PartitionConfig partition_config;
        cfg.strong(partition_config);
        partition_config.seed = seed;
        //configure strong
        internal_kaffpa_call(partition_config, suppress_output, n, vwgt, xadj, adjcwgt, adjncy, nparts, imbalance, edgecut, part);
}

void kaffpa_eco(int *n, 
                int* vwgt,
                int *xadj, 
                int* adjcwgt, 
                int *adjncy, 
                int *nparts, 
                double *imbalance,  
                bool suppress_output, 
                int seed,
                int *edgecut, 
                int *part) {

        //configure eco 
        configuration cfg;
        PartitionConfig partition_config;
        cfg.eco(partition_config);
        partition_config.seed = seed;

        internal_kaffpa_call(partition_config, suppress_output, n, vwgt, xadj, adjcwgt, adjncy, nparts, imbalance, edgecut, part);
}

void kaffpa_fast(int* n, 
                int* vwgt, 
                int* xadj, 
                int* adjcwgt, 
                int* adjncy,  
                int* nparts, 
                double* imbalance,  
                bool suppress_output, 
                int seed,
                int* edgecut, 
                int* part) {

        //configure fast 
        configuration cfg;
        PartitionConfig partition_config;
        cfg.fast(partition_config);
        partition_config.seed = seed;
        internal_kaffpa_call(partition_config, suppress_output, n, vwgt, xadj, adjcwgt, adjncy, nparts, imbalance, edgecut, part);
}


