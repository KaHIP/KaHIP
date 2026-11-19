/******************************************************************************
 * connected_blocks.cpp
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <cmath>
#include <limits>
#include <queue>
#include <vector>

#include "partition/connected_blocks.h"
#include "partition/uncoarsening/refinement/connectivity_check.h"
#include "partition/uncoarsening/refinement/mixed_refinement.h"
#include "partition/uncoarsening/refinement/quotient_graph_refinement/complete_boundary.h"

bool graph_is_connected(graph_access & G) {
        if(G.number_of_nodes() == 0) return true;
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
        return visited_count == G.number_of_nodes();
}

void enforce_connected_blocks(PartitionConfig & partition_config, graph_access & G) {
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
