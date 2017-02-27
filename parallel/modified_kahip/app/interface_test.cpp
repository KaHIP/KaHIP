/******************************************************************************
 * interface_test.cpp
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

#include <argtable2.h>
#include <iostream>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h> 

#include "data_structure/graph_access.h"
#include "graph_io.h"
#include "macros_assertions.h"
#include "parse_parameters.h"
#include "partition/coarsening/clustering/size_constraint_label_propagation.h"
#include "partition/graph_partitioner.h"
#include "partition/partition_config.h"
#include "partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"
#include "../interface/kaHIP_interface.h"


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
 
        NodeWeight largest_graph_weight = 0;
        forall_nodes(G, node) {
                largest_graph_weight += G.getNodeWeight(node);
        } endfor
        
        double epsilon = (partition_config.imbalance)/100.0;
        partition_config.upper_bound_partition      = (1+epsilon)*ceil(largest_graph_weight/(double)partition_config.k);
        partition_config.largest_graph_weight       = largest_graph_weight;
        partition_config.graph_allready_partitioned = false;
        partition_config.kway_adaptive_limits_beta  = log(largest_graph_weight);

        std::cout <<  "block weight upper bound " <<  partition_config.upper_bound_partition  << std::endl;
        
        if(partition_config.input_partition != "") {
                std::cout <<  "reading input partition" << std::endl;
                graph_io::readPartition(G, partition_config.input_partition);
                partition_config.graph_allready_partitioned  = true;
                partition_config.no_new_initial_partitioning = true;
        }

        srand(partition_config.seed);
        random_functions::setSeed(partition_config.seed);

        std::cout <<  "graph has " <<  G.number_of_nodes() <<  " nodes and " <<  G.number_of_edges() <<  " edges"  << std::endl;
        // ***************************** perform partitioning ***************************************       
        t.restart();
        graph_partitioner partitioner;
        quality_metrics qm;

        std::cout <<  "performing partitioning!"  << std::endl;

        int n = G.number_of_nodes();
        int* vwgt = G.UNSAFE_metis_style_vwgt_array();
        int* xadj = G.UNSAFE_metis_style_xadj_array();
        int* adjcwgt = G.UNSAFE_metis_style_adjwgt_array();
        int* adjncy = G.UNSAFE_metis_style_adjncy_array();
        double imbalance = 0.03;
        int* part = new int[G.number_of_nodes()];
        int edge_cut = 0;
        int nparts = partition_config.k;
        kaffpa(&n, vwgt, xadj, adjcwgt, adjncy, &nparts, &imbalance, false, 0, STRONG,& edge_cut, part);
        std::cout <<  "edge cut " <<  edge_cut << std::endl;
        forall_nodes(G, node) {
                G.setPartitionIndex(node, part[node]);
        } endfor
        
        //if(partition_config.time_limit == 0) {
                                //partitioner.perform_partitioning(partition_config, G);
        //} else {
                //PartitionID* map = new PartitionID[G.number_of_nodes()];
                //EdgeWeight best_cut = G.number_of_edges();
                //while(t.elapsed() < partition_config.time_limit) {
                        //partition_config.graph_allready_partitioned = false;
                        //partitioner.perform_partitioning(partition_config, G);
                        //EdgeWeight cut = qm.edge_cut(G);
                        //if(cut < best_cut) {
                                //best_cut = cut;
                                //forall_nodes(G, node) {
                                        //map[node] = G.getPartitionIndex(node);
                                //} endfor
                        //}
                //}

                //forall_nodes(G, node) {
                        //G.setPartitionIndex(node, map[node]);
                //} endfor

        //}

        //if( partition_config.kaffpa_perfectly_balance ) {
                //double epsilon                         = partition_config.imbalance/100.0;
                //partition_config.upper_bound_partition = (1+epsilon)*ceil(largest_graph_weight/(double)partition_config.k);

                //complete_boundary boundary(&G);
                //boundary.build();

                //cycle_refinement cr;
                //cr.perform_refinement(partition_config, G, boundary);
        //}
        //// ******************************* done partitioning *****************************************       
        //ofs.close();
        //std::cout.rdbuf(backup);
        //std::cout <<  "time spent for partitioning " << t.elapsed()  << std::endl;
       
        //// output some information about the partition that we have computed 
        //std::cout << "cut \t\t"         << qm.edge_cut(G)                 << std::endl;
        std::cout << "finalobjective  " << qm.edge_cut(G)                 << std::endl;
        //std::cout << "bnd \t\t"         << qm.boundary_nodes(G)           << std::endl;
        std::cout << "balance \t"       << qm.balance(G)                  << std::endl;
        //std::cout << "finalbalance \t"  << qm.balance(G)                  << std::endl;
        //std::cout << "max_comm_vol \t"  << qm.max_communication_volume(G) << std::endl;

        //// write the partition to the disc 
        //std::string partition("tmppartition");
        //std::stringstream noparts;
        //noparts << "tmppartition" << partition_config.k;
        //graph_io::writePartition(G, noparts.str());
        
}
