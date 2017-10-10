/******************************************************************************
 * kaHIP_interface.cpp
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

#include <iostream>
#include "kaHIP_interface.h"
#include "../lib/data_structure/graph_access.h"
#include "../lib/io/graph_io.h"
#include "../lib/tools/timer.h"
#include "../lib/tools/quality_metrics.h"
#include "../lib/tools/macros_assertions.h"
#include "../lib/tools/random_functions.h"
//#include "../lib/parallel_mh/parallel_mh_async.h"
#include "../lib/partition/uncoarsening/separator/area_bfs.h"
#include "../lib/partition/partition_config.h"
#include "../lib/partition/graph_partitioner.h"
#include "../lib/partition/uncoarsening/separator/vertex_separator_algorithm.h"
#include "../app/configuration.h"
#include "../app/balance_configuration.h"

using namespace std;

void internal_build_graph
(
    PartitionConfig & partition_config,
    const int n,
    const int* vwgt,
    const int* xadj,
    const int* adjcwgt,
    const int* adjncy,
    graph_access & G
) {
        G.build_from_metis(n, xadj, adjncy);
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

        balance_configuration bc;
        bc.configurate_balance( partition_config, G);
}

void internal_kaffpa_call
(
    PartitionConfig & partition_config,
    bool suppress_output,
    const int n,
    const int* vwgt,
    const int* xadj,
    const int* adjcwgt,
    const int* adjncy,
    const double imbalance,
    int* edgecut,
    int* part
) {

        streambuf* backup = cout.rdbuf();
        ofstream ofs;
        ofs.open("/dev/null");
        if(suppress_output) {
               cout.rdbuf(ofs.rdbuf()); 
        }

        partition_config.imbalance = 100*(imbalance);
        graph_access G;     
        internal_build_graph( partition_config, n, vwgt, xadj, adjcwgt, adjncy, G);

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

void kaffpa
(
    const int* n,
    const int* vwgt,
    const int* xadj,
    const int* adjcwgt,
    const int* adjncy,
    const int* nparts,
    const double* imbalance,
    bool suppress_output,
    int seed,
    int mode,
    int* edgecut,
    int* part
)
{
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
        internal_kaffpa_call
        (
            partition_config,
            suppress_output,
            *n,
            vwgt,
            xadj,
            adjcwgt,
            adjncy,
            *imbalance,
            edgecut,
            part
        );
}

void kaffpa_balance_NE
(
    const int* n,
    const int* vwgt,
    const int* xadj,
    const int* adjcwgt,
    const int* adjncy,
    const int* nparts,
    const double* imbalance,
    bool suppress_output,
    int seed,
    int mode,
    int* edgecut,
    int* part
)
{
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
        partition_config.balance_edges = true;
        internal_kaffpa_call
        (
            partition_config,
            suppress_output,
            *n,
            vwgt,
            xadj,
            adjcwgt,
            adjncy,
            *imbalance,
            edgecut,
            part
        );
}

void internal_nodeseparator_call
(
    PartitionConfig & partition_config,
    bool suppress_output,
    const int n,
    const int* vwgt,
    const int* xadj,
    const int* adjcwgt,
    const int* adjncy,
    const double imbalance,
    int mode,
    int* num_nodeseparator_vertices,
    int** separator
)
{

        //first perform std partitioning using KaFFPa
        streambuf* backup = cout.rdbuf();
        ofstream ofs;
        ofs.open("/dev/null");
        if(suppress_output) {
               cout.rdbuf(ofs.rdbuf()); 
        }

        // partition_config.k : already set by the caller
        partition_config.imbalance = 100*(imbalance);
        graph_access G;     
        internal_build_graph( partition_config, n, vwgt, xadj, adjcwgt, adjncy, G);
        graph_partitioner partitioner;

        area_bfs::m_deepth.resize(G.number_of_nodes());
        forall_nodes(G, node) {
                area_bfs::m_deepth[node] = 0;
        } endfor

        if( partition_config.k > 2 ) {
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
        } else {
                
                configuration cfg;
                switch( mode ) {
                        case FAST: 
                                cfg.fast_separator(partition_config);
                                break;
                        case ECO: 
                                cfg.eco_separator(partition_config);
                                break;
                        case STRONG: 
                                cfg.strong_separator(partition_config);
                                break;
                        case FASTSOCIAL: 
                                cfg.fast_separator(partition_config);
                                //cfg.fastsocial_separator(partition_config);
                                break;
                        case ECOSOCIAL: 
                                cfg.eco_separator(partition_config);
                                //cfg.ecosocial_separator(partition_config);
                                break;
                        case STRONGSOCIAL: 
                                //cfg.strongsocial_separator(partition_config);
                                cfg.strong_separator(partition_config);
                                break;
                        default: 
                                cfg.strong_separator(partition_config);
                                break;
                }       
                partition_config.mode_node_separators = true;
                partitioner.perform_partitioning(partition_config, G);
                NodeWeight ns_size = 0;
                forall_nodes(G, node) {
                        if(G.getPartitionIndex(node) == G.getSeparatorBlock()) {
                                ns_size++;
                        }
                } endfor
                *num_nodeseparator_vertices = ns_size;
                *separator = new int[*num_nodeseparator_vertices];
                unsigned int i = 0;
                forall_nodes(G, node) {
                        if(G.getPartitionIndex(node) == G.getSeparatorBlock()) {
                                (*separator)[i] = node;
                                i++;
                        }
                } endfor
        }

        ofs.close();
        cout.rdbuf(backup);
}


void node_separator
(
    const int* n,
    const int* vwgt,
    const int* xadj,
    const int* adjcwgt,
    const int* adjncy,
    const int* nparts,
    const double* imbalance,
    bool suppress_output,
    int seed,
    int mode,
    int* num_separator_vertices,
    int** separator
)
{
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

        internal_nodeseparator_call
        (
            partition_config,
            suppress_output,
            *n,
            vwgt,
            xadj,
            adjcwgt,
            adjncy,
            *imbalance,
            mode,
            num_separator_vertices,
            separator
        );
}


