/*
 * Author: Wolfgang Ost
 */

#include <algorithm>
#include <cassert>
#include <numeric>
#include <queue>
#include <unordered_map>
#include <vector>

#include "data_structure/graph_access.h"
#include "data_structure/union_find.h"
#include "node_ordering/ordering_tools.h"

void print_ordering(std::ostream &out, const std::vector<NodeID> &labels) {
        out << labels.size() << std::endl;
        for (size_t i = 0; i < labels.size(); ++i) {
                out << (i + 1) << "\t" << (labels[i] + 1) << std::endl;
        }
}

NodeWeight compute_reachable_set_size(graph_access &graph, NodeID node) {
        NodeWeight result = graph.getNodeWeight(node) - 1;
        forall_out_edges(graph, edge, node) {
                result += graph.getNodeWeight(graph.getEdgeTarget(edge));
        } endfor
        result -= graph.get_contraction_offset(node);
        return result;
}

// Compute the number of fill-edges using the algorithm in
// Rose et al. "Algorithmic Aspects of Vertex Elimination on Graphs", SIAM J. Comput., Vol. 5, No. 2, 1976 (RTL76)
Count compute_fill(graph_access &graph, const std::vector<NodeID> &ordering) {
        // Find the order of nodes by sorting node ids based on their order
        std::vector<NodeID> node_order(graph.number_of_nodes());
        std::iota(node_order.begin(), node_order.end(), 0);
        std::sort(node_order.begin(), node_order.end(),
                        [&ordering](NodeID a, NodeID b) { return ordering[a] < ordering[b]; });

        // Build a copy of 'graph' that we can modify, using monotone adjacencies
        // The monotone adjacency of x are those neighbors of x which are ordered after x.
        std::vector<std::vector<NodeID>> adjacencies;
        Count edge_count = 0;
        forall_nodes(graph, node) {
                adjacencies.push_back({});
                forall_out_edges(graph, edge, node) {
                        auto target = graph.getEdgeTarget(edge);
                        if (ordering[target] > ordering[node]) {
                                adjacencies.back().push_back(target);
                                edge_count++;
                        }
                } endfor
        } endfor

        std::vector<bool> test(graph.number_of_nodes(), false);
        for (const auto node: node_order) {             // node = i in RTL76
                NodeID k = graph.number_of_nodes();            // smallest ordering in the neighborhood of node that's greater than node
                // eliminate duplicates in A(v) and compute m(v), v = node
                // A(v): monotone adjacency of v
                // m(v): node ordered at position k
                for (auto it = adjacencies[node].begin(); it != adjacencies[node].end(); ) {
                        auto neighbor = *it;            // neighbor = w in RTL76
                        if (test[ordering[neighbor]]) {
                                it = adjacencies[node].erase(it);
                        } else {
                                test[ordering[neighbor]] = true;
                                k = std::min(k, ordering[neighbor]);
                                ++it;
                        }
                }
                auto m = node_order[k];         // the node ordered in the k-th place
                for (auto it = adjacencies[node].begin(); it != adjacencies[node].end(); ++it) {
                        auto neighbor = *it;
                        test[ordering[neighbor]] = false;
                        if (neighbor != m) {
                                adjacencies[m].push_back(neighbor);
                        }
                }
        }
        Count fill_edge_count = 0;
        for (const auto &adj: adjacencies) {
                fill_edge_count += adj.size();
        }
        return fill_edge_count - edge_count;
}
