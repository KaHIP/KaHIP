/******************************************************************************
 * kaffpa.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <argtable3.h>
#include <iostream>
#include <math.h>
#include <queue>
#include <set>
#ifndef _WIN32
#include <regex.h>
#endif
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
#ifndef _WIN32
#include "mmap_graph_io.h"
#endif
#include "parse_parameters.h"
#include "partition/graph_partitioner.h"
#include "partition/partition_config.h"
#include "partition/uncoarsening/refinement/connectivity_check.h"
#include "partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.h"
#include "partition/uncoarsening/refinement/mixed_refinement.h"
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
#ifdef _WIN32
        ofs.open("NUL");
#else
        ofs.open("/dev/null");
#endif
        if(suppress_output) {
                std::cout.rdbuf(ofs.rdbuf()); 
        }

        partition_config.LogDump(stdout);
        graph_access G;     

        timer t;
#ifndef _WIN32
        if (partition_config.use_mmap_io) {
                kahip::mmap_io::graph_from_metis_file(G, graph_filename);
        } else
#endif
        {
                graph_io::readGraphWeighted(G, graph_filename);
        }
        std::cout << "io time: " << t.elapsed()  << std::endl;

        G.set_partition_count(partition_config.k); 

        balance_configuration bc;
        bc.configurate_balance( partition_config, G);

        srand(partition_config.seed);
        random_functions::setSeed(partition_config.seed);

        std::cout <<  "graph has " <<  G.number_of_nodes() <<  " nodes and " <<  G.number_of_edges() <<  " edges"  << std::endl;

        if(partition_config.connected_blocks) {
                std::vector<bool> visited(G.number_of_nodes(), false);
                std::queue<NodeID> bfs_queue;
                visited[0] = true;
                bfs_queue.push(0);
                NodeID visited_count = 1;
                while(!bfs_queue.empty()) {
                        NodeID v = bfs_queue.front(); bfs_queue.pop();
                        forall_out_edges(G, e, v) {
                                NodeID u = G.getEdgeTarget(e);
                                if(!visited[u]) { visited[u] = true; visited_count++; bfs_queue.push(u); }
                        } endfor
                }
                if(visited_count < G.number_of_nodes()) {
                        std::cout << "WARNING: input graph is disconnected, connected blocks cannot be guaranteed." << std::endl;
                }
        }

        // ***************************** perform partitioning ***************************************
        t.restart();
        graph_partitioner partitioner;
        quality_metrics qm;

        bool do_connected_blocks = partition_config.connected_blocks;

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
        partition_config.connected_blocks = do_connected_blocks;
        if(do_connected_blocks) {
                double epsilon = partition_config.imbalance/100.0;
                partition_config.upper_bound_partition = (1+epsilon)*ceil(partition_config.largest_graph_weight/(double)partition_config.k);

                // Iterate: eliminate disconnected components, then rebalance
                for(int iter = 0; iter < 5; iter++) {
                        // Step 1: Eliminate disconnected components (like METIS EliminateComponents)
                        bool had_disconnected = false;
                        for(PartitionID block = 0; block < partition_config.k; block++) {
                                std::vector<int> comp_id(G.number_of_nodes(), -1);
                                int num_comps = 0;
                                std::vector<NodeWeight> comp_weight;
                                std::vector<std::vector<NodeID>> comp_nodes;

                                forall_nodes(G, node) {
                                        if(G.getPartitionIndex(node) == block && comp_id[node] == -1) {
                                                NodeWeight cw = 0;
                                                std::vector<NodeID> nodes;
                                                std::queue<NodeID> q;
                                                q.push(node); comp_id[node] = num_comps;
                                                while(!q.empty()) {
                                                        NodeID v = q.front(); q.pop();
                                                        cw += G.getNodeWeight(v);
                                                        nodes.push_back(v);
                                                        forall_out_edges(G, e, v) {
                                                                NodeID u = G.getEdgeTarget(e);
                                                                if(G.getPartitionIndex(u) == block && comp_id[u] == -1) {
                                                                        comp_id[u] = num_comps; q.push(u);
                                                                }
                                                        } endfor
                                                }
                                                comp_weight.push_back(cw);
                                                comp_nodes.push_back(nodes);
                                                num_comps++;
                                        }
                                } endfor

                                if(num_comps <= 1) continue;
                                had_disconnected = true;

                                int largest = 0;
                                for(int c = 1; c < num_comps; c++) {
                                        if(comp_weight[c] > comp_weight[largest]) largest = c;
                                }

                                // Move small components to lightest adjacent block
                                for(int c = 0; c < num_comps; c++) {
                                        if(c == largest) continue;
                                        std::vector<NodeWeight> bw(partition_config.k, 0);
                                        forall_nodes(G, n) { bw[G.getPartitionIndex(n)] += G.getNodeWeight(n); } endfor

                                        PartitionID target = block;
                                        NodeWeight min_w = std::numeric_limits<NodeWeight>::max();
                                        for(NodeID n : comp_nodes[c]) {
                                                forall_out_edges(G, e, n) {
                                                        PartitionID ab = G.getPartitionIndex(G.getEdgeTarget(e));
                                                        if(ab != block && bw[ab] < min_w) { min_w = bw[ab]; target = ab; }
                                                } endfor
                                        }
                                        if(target != block) {
                                                for(NodeID n : comp_nodes[c]) G.setPartitionIndex(n, target);
                                        }
                                }
                        }

                        if(!had_disconnected && iter > 0) break;

                        // Step 2: Greedy rebalance - move non-articulation boundary
                        // nodes from overweight blocks to lightest underweight neighbor
                        bool rebal_progress = true;
                        while(rebal_progress) {
                                rebal_progress = false;
                                std::vector<NodeWeight> bw(partition_config.k, 0);
                                forall_nodes(G, n) { bw[G.getPartitionIndex(n)] += G.getNodeWeight(n); } endfor

                                forall_nodes(G, node) {
                                        PartitionID from = G.getPartitionIndex(node);
                                        if(bw[from] <= partition_config.upper_bound_partition) continue;

                                        // Find lightest underweight adjacent block
                                        PartitionID target = from;
                                        NodeWeight min_w = std::numeric_limits<NodeWeight>::max();
                                        forall_out_edges(G, e, node) {
                                                PartitionID ab = G.getPartitionIndex(G.getEdgeTarget(e));
                                                if(ab != from && bw[ab] + G.getNodeWeight(node) <= partition_config.upper_bound_partition && bw[ab] < min_w) {
                                                        min_w = bw[ab];
                                                        target = ab;
                                                }
                                        } endfor

                                        if(target == from) continue;
                                        if(would_disconnect_block(G, node, from)) continue;

                                        bw[from] -= G.getNodeWeight(node);
                                        bw[target] += G.getNodeWeight(node);
                                        G.setPartitionIndex(node, target);
                                        rebal_progress = true;
                                } endfor
                        }
                }

                // Final refinement pass with connectivity guard active
                partition_config.connected_blocks = true;
                complete_boundary final_boundary(&G);
                final_boundary.build();
                refinement* final_refine = new mixed_refinement();
                final_refine->perform_refinement(partition_config, G, final_boundary);
                delete final_refine;
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
#ifndef NDEBUG
        if(partition_config.connected_blocks) {
                for(PartitionID block = 0; block < partition_config.k; block++) {
                        NodeID start = std::numeric_limits<NodeID>::max();
                        NodeID block_size = 0;
                        forall_nodes(G, node) {
                                if(G.getPartitionIndex(node) == block) {
                                        if(start == std::numeric_limits<NodeID>::max()) start = node;
                                        block_size++;
                                }
                        } endfor

                        if(block_size == 0) continue;

                        std::vector<bool> visited(G.number_of_nodes(), false);
                        std::queue<NodeID> q;
                        visited[start] = true;
                        q.push(start);
                        NodeID reached = 1;
                        while(!q.empty()) {
                                NodeID v = q.front(); q.pop();
                                forall_out_edges(G, e, v) {
                                        NodeID u = G.getEdgeTarget(e);
                                        if(!visited[u] && G.getPartitionIndex(u) == block) {
                                                visited[u] = true;
                                                reached++;
                                                q.push(u);
                                        }
                                } endfor
                        }
                        ASSERT_EQ(reached, block_size);
                }
        }
#endif
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
