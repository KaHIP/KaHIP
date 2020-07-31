/*
 * Author: Wolfgang Ost
 */

#include <algorithm>
#include <functional>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

#include "node_ordering/min_degree_ordering.h"
#include "node_ordering/ordering_tools.h"
#include "tools/macros_assertions.h"

MinDegree::MinDegree(graph_access * const graph, const std::vector<NodeID> &halo_nodes) : graph(graph),
                                                   memberships(graph->number_of_nodes()),
                                                   node_links(graph->number_of_nodes()),
                                                   link_weights(graph->number_of_nodes()),
                                                   is_halo_node(graph->number_of_nodes(), false) {
        for (const auto node: halo_nodes) {
                is_halo_node[node] = true;
        }
        std::iota(node_links.begin(), node_links.end(), 0);
        forall_nodes((*graph), node) {
                link_weights[node] = graph->getNodeWeight(node);
        } endfor
}

void MinDegree::initialize_cliques() {
        forall_nodes((*graph), node) {
                if (!is_halo_node[node]) {
                        current_nodes.insert(node);
                        // Instead of the node degree use the adjusted reachable set size
                        // We use the negative degree, because bucket_pq sorts by maximum
                        degree_pq.insert(node, -(Gain)compute_reachable_set_size(*graph, node));
                }
                // Make each edge a clique
                forall_out_edges((*graph), edge, node) {
                        auto target = graph->getEdgeTarget(edge);
                        if (!degree_pq.contains(target) || is_halo_node[target]) {
                                cliques.push_back(clique( { node, target } ));
                                memberships[node].push_back(cliques.size()-1);
                                memberships[target].push_back(cliques.size()-1);
                        }
                } endfor
        } endfor
}

void MinDegree::perform_ordering(std::vector<NodeID> &labels) {
        initialize_cliques();

        NodeID order = 0;
        // Eliminate nodes until no nodes are left
        while (!current_nodes.empty()) {
                if (degree_pq.empty()) {
                        break;
                }

                // Select the node with the smallest degree
                auto min_degree = -degree_pq.maxValue();
                auto min_degree_node = degree_pq.deleteMax();
                if (min_degree > 0 && !memberships[min_degree_node].empty()) {
                        auto merged_clique_index = eliminate_node(min_degree_node);
                        // Label nodes indistinguishable from min_degree_node before degree update
                        label_node(min_degree_node, labels, order);
                        update_node_degrees(min_degree_node, merged_clique_index);
                } else {
                        label_node(min_degree_node, labels, order);
                }

                // remove the minimum degree node and order it next
                current_nodes.erase(current_nodes.find(min_degree_node));
                // If min_degree_node was a neighbor of an eliminated node,
                // it no longer is.
                elimination_neighbors.remove(min_degree_node);
        }
        forall_nodes((*graph), node) {
                if (is_halo_node[node]) {
                        label_node(node, labels, order);
                }
        } endfor
}

size_t MinDegree::eliminate_node(NodeID node) {
        auto merged_clique_idx = memberships[node].front();
        auto& merged_clique = cliques[merged_clique_idx];
        for (size_t idx = 1; idx < memberships[node].size(); ++idx) {
                auto clique_idx = memberships[node][idx];
                merged_clique.merge(cliques[clique_idx]);
        }
        merged_clique.remove(node);
        return merged_clique_idx;
}

void MinDegree::update_node_degrees(NodeID eliminated_node, size_t merged_clique_idx) {
        auto& merged_clique = cliques[merged_clique_idx];
        // Remove all cliques containing the min degree node
        // from the clique list of each node that needs to be updated
        // and update the node degrees

        new_indistinguishable_nodes.clear();
        clique prototype;
        bool found_prototype = false;
        // Place the updated nodes in the elimination neighbors list, before actually updating the degrees
        for (const auto node: merged_clique) {
                elimination_neighbors.insert(node);
        }
        for (const auto node: merged_clique) {
                if (is_halo_node[node]) {
                        // We don't care about the correct degree of halo nodes
                        // or which cliques they are a member of, so we update them
                        continue;
                }
                // We need the set of all neighbors of 'node' to update the degree
                // 'isec' is the intersection of all cliques 'node' belongs to,
                // so we can just sum up the weights of those nodes later.
                clique isec;
                auto it = memberships[node].begin();
                while (it != memberships[node].end()) {
                        // Remove cliques that are a subset of merged_clique (element absorption)
                        // and cliques that contain the eliminated node
                        if (cliques[*it].contains(eliminated_node)) {
                                // Swap and pop for faster removal
                                *it = memberships[node].back();
                                memberships[node].pop_back();
                        } else if (merged_clique.contains_all_of(cliques[*it].begin(), cliques[*it].end())) {
                                *it = memberships[node].back();
                                memberships[node].pop_back();
                        } else {
                                isec.merge(cliques[*it]);
                                ++it;
                        }
                }
                memberships[node].push_back(merged_clique_idx);
                isec.merge(merged_clique);

                // Test if 'node' has neighbors only in the 'elimination_neighbors' set,
                // i.e., if the neighbors of 'node' have all been neighbors of an eliminated node before
                // TODO: comparing against the first isec gives indistinguishable nodes, but it's not a great solution
                if (elimination_neighbors.contains_all_of(isec.begin(), isec.end())) {
                        if (!found_prototype) {
                                prototype = isec;
                                new_indistinguishable_nodes.push_back(node);
                                found_prototype = true;
                        } else if (isec == prototype) {
                                new_indistinguishable_nodes.push_back(node);
                        }
                }

                // Perform the degree update.
                // 'isec' contains the neighbors of 'node' and 'node' itself.
                // We can ignore the contraction offset here, because
                // it becomes zero as soon as a neighbor of 'node'
                // is eliminated.
                Gain node_degree = 0;
                for (auto neighbor: isec) {
                        node_degree += link_weights[neighbor];
                }
                node_degree -= 1;
                degree_pq.changeKey(node, -node_degree);
        }
        update_indistinguishable_nodes();
}

void MinDegree::update_indistinguishable_nodes() {
        if (new_indistinguishable_nodes.size() > 1) {
                auto representative = new_indistinguishable_nodes[0];
                auto link = representative;
                for (unsigned int i = 1; i < new_indistinguishable_nodes.size(); ++i) {
                        // link from the representative to the next indistinguishable node
                        // if the representative already links to some other node, link from
                        // the first node that links to itself
                        while (node_links[link] != link) {
                                link = node_links[link];
                        }
                        link = node_links[link] = new_indistinguishable_nodes[i];
                        
                        // Remove the non-representative nodes from all relevant node sets
                        // Remove the node from the current_nodes
                        current_nodes.erase(current_nodes.find(new_indistinguishable_nodes[i]));
                        degree_pq.deleteNode(new_indistinguishable_nodes[i]);
                        // Remove the node from all cliques it's a member of
                        for (auto &clique_id: memberships[new_indistinguishable_nodes[i]]) {
                                cliques[clique_id].remove(new_indistinguishable_nodes[i]);
                        }
                        // Remove the node from elimination_neighbors?
                        elimination_neighbors.remove(new_indistinguishable_nodes[i]);
                }
                // Compute total weight of indistinguishable nodes
                link = representative;
                link_weights[representative] = graph->getNodeWeight(link);
                while (node_links[link] != link) {
                        link = node_links[link];
                        link_weights[representative] += graph->getNodeWeight(link);
                }
        }
}

void MinDegree::label_node(NodeID node, std::vector<NodeID> &labels, NodeID &order) {
        auto link = node;
        // Order the node itself
        labels[link] = order;
        order += 1;
        // Order indistinguishable nodes
        while (node_links[link] != link) {
                link = node_links[link];
                labels[link] = order;
                order += 1;
        }
}
