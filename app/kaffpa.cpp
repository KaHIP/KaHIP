/******************************************************************************
 * kaffpa.cpp 
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

#include <argtable2.h>
#include <iostream>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h> 

#include "balance_configuration.h"
#include "data_structure/graph_access.h"
#include "data_structure/matrix/normal_matrix.h"
#include "data_structure/matrix/online_distance_matrix.h"
#include "graph_io.h"
#include "macros_assertions.h"
#include "mapping/mapping_algorithms.h"
#include "parse_parameters.h"
#include "partition/graph_partitioner.h"
#include "partition/partition_config.h"
#include "partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"

int main(int argn, char **argv) {

        PartitionConfig partition_config;
        std::string graph_filename;

        bool is_graph_weighted = false;
        bool suppress_output   = false;
        bool recursive         = false;

        int ret_code = parse_parameters(argn, argv, 
                        partition_config, 
                        graph_filename, 
                        is_graph_weighted, 
                        suppress_output, recursive); 

        if(ret_code) {
                return 0;
        }

        std::streambuf* backup = std::cout.rdbuf();
        std::ofstream ofs;
        ofs.open("/dev/null");
        if(suppress_output) {
                std::cout.rdbuf(ofs.rdbuf()); 
        }

        partition_config.LogDump(stdout);
        graph_access G;     

        timer t;
        graph_io::readGraphWeighted(G, graph_filename);
        std::cout << "io time: " << t.elapsed()  << std::endl;

        G.set_partition_count(partition_config.k); 

        balance_configuration bc;
        bc.configurate_balance( partition_config, G);

        srand(partition_config.seed);
        random_functions::setSeed(partition_config.seed);

        std::cout <<  "graph has " <<  G.number_of_nodes() <<  " nodes and " <<  G.number_of_edges() <<  " edges"  << std::endl;
        // ***************************** perform partitioning ***************************************       
        t.restart();
        graph_partitioner partitioner;
        quality_metrics qm;

        std::cout <<  "performing partitioning!"  << std::endl;
        if(partition_config.time_limit == 0) {
                partitioner.perform_partitioning(partition_config, G);
        } else {
                PartitionID* map = new PartitionID[G.number_of_nodes()];
                EdgeWeight best_cut = std::numeric_limits<EdgeWeight>::max();
                while(t.elapsed() < partition_config.time_limit) {
                        partition_config.graph_allready_partitioned = false;
                        partitioner.perform_partitioning(partition_config, G);
                        EdgeWeight cut = qm.edge_cut(G);
                        if(cut < best_cut) {
                                best_cut = cut;
                                forall_nodes(G, node) {
                                        map[node] = G.getPartitionIndex(node);
                                } endfor
                        }
                }

                forall_nodes(G, node) {
                        G.setPartitionIndex(node, map[node]);
                } endfor
        }

        if( partition_config.kaffpa_perfectly_balance ) {
                double epsilon                         = partition_config.imbalance/100.0;
                partition_config.upper_bound_partition = (1+epsilon)*ceil(partition_config.largest_graph_weight/(double)partition_config.k);

                complete_boundary boundary(&G);
                boundary.build();

                cycle_refinement cr;
                cr.perform_refinement(partition_config, G, boundary);
        }
        ofs.close();
        std::cout.rdbuf(backup);
        std::cout <<  "time spent for partitioning " << t.elapsed()  << std::endl;

        int qap = 0;
        if(partition_config.enable_mapping) {
                std::cout <<  "performing mapping!"  << std::endl;
                //check if k is a power of 2 
                bool power_of_two = (partition_config.k & (partition_config.k-1)) == 0;
                std::vector< NodeID > perm_rank(partition_config.k);
                graph_access C;
                complete_boundary boundary(&G);
                boundary.build();
                boundary.getUnderlyingQuotientGraph(C);

                forall_nodes(C, node) {
                        C.setNodeWeight(node, 1);
                } endfor

                if(!power_of_two ) {
                        t.restart();
                        mapping_algorithms ma;
                        if( partition_config.distance_construction_algorithm != DIST_CONST_HIERARCHY_ONLINE) {
                                normal_matrix D(partition_config.k, partition_config.k);
                                ma.construct_a_mapping(partition_config, C, D, perm_rank);
                                std::cout <<  "time spent for mapping " << t.elapsed()  << std::endl;
                                qap = qm.total_qap(C, D, perm_rank );
                        } else {
                                online_distance_matrix D(partition_config.k, partition_config.k);
                                D.setPartitionConfig(partition_config);
                                ma.construct_a_mapping(partition_config, C, D, perm_rank);
                                std::cout <<  "time spent for mapping " << t.elapsed()  << std::endl;
                                qap = qm.total_qap(C, D, perm_rank );
                        }
                } else {
                        std::cout <<  "number of processors is a power of two, so no mapping algorithm is performed (identity is best)"  << std::endl;
                        std::cout <<  "time spent for mapping " << 0 << std::endl;
                        for( unsigned i = 0; i < perm_rank.size(); i++) {
                                perm_rank[i] = i;
                        }

                        online_distance_matrix D(partition_config.k, partition_config.k);
                        D.setPartitionConfig(partition_config);
                        qap = qm.total_qap(C, D, perm_rank );
                }

                // solution check 
                std::vector< NodeID > tbsorted = perm_rank;
                std::sort( tbsorted.begin(), tbsorted.end() );
                for( unsigned int i = 0; i < tbsorted.size(); i++) {
                        if( tbsorted[i] != i ) {
                                std::cout <<  "solution is NOT a permutation. Please report this."  << std::endl;
                                std::cout <<  tbsorted[i] <<  " " << i   << std::endl;
                                exit(0);
                        }
                }

                forall_nodes(G, node) {
                        G.setPartitionIndex(node, perm_rank[G.getPartitionIndex(node)]);
                } endfor
        }
        // ******************************* done partitioning *****************************************       
        // output some information about the partition that we have computed 
        std::cout << "cut \t\t"         << qm.edge_cut(G)                 << std::endl;
        std::cout << "finalobjective  " << qm.edge_cut(G)                 << std::endl;
        std::cout << "bnd \t\t"         << qm.boundary_nodes(G)           << std::endl;
        std::cout << "balance \t"       << qm.balance(G)                  << std::endl;
        std::cout << "max_comm_vol \t"  << qm.max_communication_volume(G) << std::endl;
        if(partition_config.enable_mapping) std::cout <<  "quadratic assignment objective J(C,D,Pi') = " << qap << std::endl;

        // write the partition to the disc 
        std::stringstream filename;
        if(!partition_config.filename_output.compare("")) {
                filename << "tmppartition" << partition_config.k;
        } else {
                filename << partition_config.filename_output;
        }

        graph_io::writePartition(G, filename.str());

}
