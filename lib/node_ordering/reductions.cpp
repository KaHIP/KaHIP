/*
 * Author: Wolfgang Ost
 */

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <queue>
#include <unordered_set>
#include <utility>

#include "data_structure/priority_queues/bucket_pq.h"
#include "node_ordering/ordering_tools.h"
#include "node_ordering/reductions.h"
#include "tools/macros_assertions.h"
#include "tools/timer.h"

#include "io/graph_io.h"

/******************************/
/* NODE CONTRACTION FUNCTIONS */
/******************************/

// Input:
//  - 'graph_before':           the graph before contraction
//  - 'node_groups':            the groups of nodes to be contracted
//
// Outputs:
//  - 'graph_after':            the graph after contraction
//  - 'mapping':                mapping from nodes in 'graph_after' to nodes in 'graph_before'
void contract_nodes(graph_access &graph_before, graph_access &graph_after,
                    const std::vector<std::vector<NodeID>> &node_groups,
                    std::unordered_map<NodeID, std::vector<NodeID>> &mapping) {
        graph_after.start_construction(node_groups.size(), graph_before.number_of_edges());

        std::vector<NodeID> reverse_map(graph_before.number_of_nodes(), 0);
        for (const auto &group: node_groups) {
                auto new_node_id = graph_after.new_node();
                NodeWeight weight = 0;
                NodeWeight offset = 0;
                for (auto old_node_id: group) {
                        reverse_map[old_node_id] = new_node_id;
                        weight += graph_before.getNodeWeight(old_node_id);
                        offset += graph_before.get_contraction_offset(old_node_id);
                }
                graph_after.setNodeWeight(new_node_id, weight);
                graph_after.set_contraction_offset(new_node_id, offset);
                mapping.insert({new_node_id, group});
        }

        std::vector<NodeID> visited_by(graph_before.number_of_nodes(), -1); // assuming that there are less than 2^32-1 nodes
        std::vector<EdgeID> target_edges(graph_before.number_of_nodes(), -1); // new_edge_ids for a target
        for (NodeID contracted_node_id = 0; contracted_node_id < node_groups.size(); ++contracted_node_id) {
                const auto &group = node_groups[contracted_node_id];
                // All out-edges from each node in the group are potential edges in the reduced graph
                for (auto node: group) {
                        forall_out_edges(graph_before, edge, node) {
                                auto target = graph_before.getEdgeTarget(edge);
                                // don't add edges between contracted nodes
                                if (reverse_map[target] != contracted_node_id) {
                                        if (visited_by[reverse_map[target]] != contracted_node_id) {
                                                // this edge hasn't been added yet
                                                auto new_edge_id = graph_after.new_edge(contracted_node_id, reverse_map[target]);
                                                target_edges[reverse_map[target]] = new_edge_id;
                                                graph_after.setEdgeWeight(new_edge_id, graph_before.getEdgeWeight(edge));
                                                visited_by[reverse_map[target]] = contracted_node_id;
                                        } else {
                                                // add to the edge weight of an existing edge
                                                auto edge_id = target_edges[reverse_map[target]];
                                                auto current_edge_weight = graph_after.getEdgeWeight(edge_id);
                                                graph_after.setEdgeWeight(edge_id, current_edge_weight + graph_before.getEdgeWeight(edge));
                                        }
                                }
                        } endfor
                }
        }

        graph_after.finish_construction();
}

/**********************************************/
/* SIMPLICIAL NODES/ISOLATED CLIQUE REDUCTION */
/**********************************************/

// Helper function to test if a node and its neighborhood form a clique
// G:                   input graph
// node:                node to be tested
// actual_degree:       degree of 'node', ignoring removed nodes
// labels:              vector to count neighbors
// removed:             mark simplicial nodes as removed
bool clique_test(graph_access &G, NodeID node, Gain actual_degree, std::vector<short> &labels, const std::vector<bool> &removed) {
        // To test if 'node' is simplicial, we iterate through its neighbors.
        // In each iteration, we
        //  - mark the currently visited node
        //  - count the neighbors with a mark
        //  - if the count equals the number of already visited neighbors,
        //      then the node could be simplicial
        //  - otherwise, there is a pair of neighbors that's not connected by an edge

        // Update labels
        bool simplicial = true;
        Count count;
        // Neighbors visited so far
        Count max_count = 0;
        forall_out_edges(G, edge, node) {
                auto neighbor = G.getEdgeTarget(edge);
                if (!removed[neighbor]) {
                        labels[neighbor] = 1;
                        if (G.getNodeDegree(neighbor) < actual_degree ||
                            G.get_contraction_offset(neighbor) != 0) {
                                simplicial = false;
                                break;
                        }
                        count = 0;
                        forall_out_edges(G, edge2, neighbor) {
                                count += labels[G.getEdgeTarget(edge2)];
                                if (count == max_count) {
                                        max_count++;
                                        // 'node' might still be simplicial,
                                        // so continue the outer loop
                                        goto still_simplicial;
                                }
                        } endfor
                        // 'neighbor' is not adjacent to 'max_count' previously visited
                        // neighbors of 'node', and thus, 'node' is not simplicial.
                        simplicial = false;
                        break;
                }
                still_simplicial:;
        } endfor
        
        // Reset labels
        forall_out_edges(G, edge, node) {
                labels[G.getEdgeTarget(edge)] = 0;
        } endfor

        return simplicial;
}

class bucket_sorter {
public:
        bucket_sorter(graph_access &graph,
                      Count degree_limit) : buckets(graph.getMaxDegree() + 1, std::vector<NodeID>{}),
                                            locations(graph.number_of_nodes(), {-1, 0}),
                                            min_ptr(graph.getMaxDegree() + 1),
                                            num_elements(0) {
                forall_nodes(graph, node) {
                        Count degree = graph.getNodeDegree(node);
                        if (degree <= degree_limit) {
                                locations[node] = {degree, buckets[degree].size()};
                                buckets[degree].push_back(node);
                                min_ptr = std::min(min_ptr, degree);
                                num_elements++;
                        }
                } endfor
        }

        bool contains(NodeID node) const {
                return locations[node].first >= 0;
        }

        void decreaseKey(NodeID node, Count decrease) {
                remove_from_bucket(locations[node].first, locations[node].second);
                locations[node].first -= decrease;
                locations[node].second = buckets[locations[node].first].size();
                buckets[locations[node].first].push_back(node);
                if (locations[node].first < (int)min_ptr) {
                        min_ptr = locations[node].first;
                }
        }

        NodeID deleteMin() {
                auto min_node = buckets[min_ptr].back();
                buckets[min_ptr].pop_back();
                locations[min_node].first = -1;
                num_elements--;
                if (num_elements > 0) {
                        while (buckets[min_ptr].empty()) {
                                min_ptr++;
                        }
                } else {
                        min_ptr = buckets.size() + 1;
                }
                return min_node;
        }
        
        Count minValue() const {
                return min_ptr;
        }

        Count size() const {
                return num_elements;
        }

private:
        std::vector<std::vector<NodeID>> buckets;
        // bucket and index in bucket
        std::vector<std::pair<int, size_t>> locations;
        Count min_ptr;
        Count num_elements;

        void remove_from_bucket(int bucket, size_t index) {
                locations[buckets[bucket].back()].second = index;
                std::swap(buckets[bucket][index], buckets[bucket].back());
                buckets[bucket].pop_back();
        }
};

void SimplicialNodeReduction::apply() {
        bucket_sorter degree_queue(graph_before, degree_limit);

        graph_after.start_construction(graph_before.number_of_nodes(), graph_before.number_of_edges());

        // mapping from nodes of graph_before to nodes of graph_after
        std::vector<NodeID> reverse_mapping(graph_before.number_of_nodes(), 0);
        std::vector<bool> remove(graph_before.number_of_nodes(), false);
        std::vector<short> test_labels(graph_before.number_of_nodes(), 0);
        while (degree_queue.size() > 0) {
                auto current_degree = degree_queue.minValue();
                NodeID node_to_test = degree_queue.deleteMin();
                if (current_degree <= 1 || clique_test(graph_before, node_to_test, current_degree, test_labels, remove)) {
                        remove[node_to_test] = true;
                        label_first.push_back(node_to_test);

                        // Update the degree of all neighbors of the eliminated nodes.
                        // The eliminated nodes are all only a member of one clique,
                        // so the neighbors of node_to_test are all we need to update.
                        forall_out_edges(graph_before, edge, node_to_test) {
                                auto target = graph_before.getEdgeTarget(edge);
                                if (degree_queue.contains(target)) {
                                        degree_queue.decreaseKey(target, 1);
                                }
                        } endfor
                }
        }

        forall_nodes(graph_before, node) {
                if (!remove[node]) {
                        auto new_node_id = graph_after.new_node();
                        graph_after.setNodeWeight(new_node_id, graph_before.getNodeWeight(node));
                        graph_after.set_contraction_offset(new_node_id, graph_before.get_contraction_offset(node));
                        mapping.push_back(node);
                        reverse_mapping[node] = new_node_id;
                }
        } endfor

        // Copy edges
        for (NodeID new_node_id = 0; new_node_id < graph_before.number_of_nodes() - label_first.size(); ++new_node_id) {
                forall_out_edges(graph_before, edge, mapping[new_node_id]) {
                        auto target = graph_before.getEdgeTarget(edge);
                        if (!remove[target]) {
                                auto new_edge_id = graph_after.new_edge(new_node_id, reverse_mapping[target]);
                                graph_after.setEdgeWeight(new_edge_id, graph_before.getEdgeWeight(edge));
                        }
                } endfor
        }

        graph_after.finish_construction();
}

void SimplicialNodeReduction::map(std::vector<NodeID> &reduced_label, std::vector<NodeID> &new_label) const {
        new_label.resize(graph_before.number_of_nodes());
        NodeID order = 0;
        for (auto node: label_first) {
                new_label[node] = order;
                order++;
        }
        for (size_t i = 0; i < reduced_label.size(); ++i) {
                new_label[mapping[i]] = reduced_label[i] + order;
        }
}

/***************************/
/* INDISTINGUISHABLE NODES */
/***************************/

size_t open_neighborhood_hash(graph_access &graph, NodeID node) {
        size_t result = 0;
        forall_out_edges(graph, edge, node) {
                result += graph.getEdgeTarget(edge);
        } endfor
        return result;
}

size_t closed_neighborhood_hash(graph_access &graph, NodeID node) {
        return open_neighborhood_hash(graph, node) + node;
}

void IndistinguishableNodeReduction::apply() {
        // Vector of hashes of closed neighborhoods
        std::vector<size_t> hash_vector;
        hash_vector.reserve(graph_before.number_of_nodes());
        forall_nodes(graph_before, node) {
                hash_vector.push_back(closed_neighborhood_hash(graph_before, node));
        } endfor

        // Sets of indistinguishable nodes
        std::vector<std::vector<NodeID>> node_groups;
        node_groups.reserve(graph_before.number_of_nodes());

        // Label nodes to compare neighborhoods
        std::vector<short> labels(graph_before.number_of_nodes(), 0);

        // Set to true, if a node has been contracted
        std::vector<bool> contracted(graph_before.number_of_nodes(), false);

        // Store adjusted degrees for each node
        std::vector<NodeWeight> degrees(graph_before.number_of_nodes(), 0);
        forall_nodes(graph_before, node) {
                degrees[node] = compute_reachable_set_size(graph_before, node);
        } endfor

        forall_nodes(graph_before, node_a) {
                if (contracted[node_a]) {
                        continue;
                }
                node_groups.push_back({node_a});
                contracted[node_a] = true;

                forall_out_edges(graph_before, edge, node_a) {
                        labels[graph_before.getEdgeTarget(edge)] = 1;
                } endfor
                labels[node_a] = 1;
                
                forall_out_edges(graph_before, edge, node_a) {
                        auto node_b = graph_before.getEdgeTarget(edge);
                        // node_b hast been visited before
                        // or the hashes are not equal
                        // or the actual degrees are not equal
                        // or the adjusted degrees are not equal
                        if (contracted[node_b] ||
                            hash_vector[node_a] != hash_vector[node_b] ||
                            graph_before.getNodeDegree(node_a) != graph_before.getNodeDegree(node_b) ||
                            degrees[node_a] != degrees[node_b]) {
                                continue;
                        }
                        
                        EdgeWeight count = 0;
                        forall_out_edges(graph_before, edge_b, node_b) {
                                count += labels[graph_before.getEdgeTarget(edge_b)];
                        } endfor
                        if (count == graph_before.getNodeDegree(node_a) && count == graph_before.getNodeDegree(node_b)) {
                                node_groups.back().push_back(node_b);
                                contracted[node_b] = true;
                        }
                } endfor

                forall_out_edges(graph_before, edge, node_a) {
                        labels[graph_before.getEdgeTarget(edge)] = 0;
                } endfor
                labels[node_a] = 0;
        } endfor

        contract_nodes(graph_before, graph_after, node_groups, mapping);

        // Match reachable set size in reduced and original graph
        for (NodeID new_node_id = 0; new_node_id < node_groups.size(); ++new_node_id) {
                if (graph_before.get_contraction_offset(node_groups[new_node_id][0]) == 0) {
                        continue;
                }
                auto reach = compute_reachable_set_size(graph_before, node_groups[new_node_id][0]);
                auto new_reach_without_offset = graph_after.getNodeWeight(new_node_id) - 1;
                forall_out_edges(graph_after, edge, new_node_id) {
                        new_reach_without_offset += graph_after.getNodeWeight(graph_after.getEdgeTarget(edge));
                } endfor
                auto offset = new_reach_without_offset - reach;
                graph_after.set_contraction_offset(new_node_id, offset);
        }
}

void IndistinguishableNodeReduction::map(std::vector<NodeID> &reduced_label, std::vector<NodeID> &new_label) const {
        new_label.resize(graph_before.number_of_nodes());
        std::vector<NodeID> orderings(reduced_label.size(), 0);
        for (size_t node = 0; node < reduced_label.size(); node++) {
                orderings[reduced_label[node]] = node;
        }
        NodeID offset = 0;
        for (size_t order = 0; order < orderings.size(); ++order) {
                auto node = orderings[order];
                NodeID count = 0;
                for (const auto original_node: mapping.at(node)) {
                        new_label[original_node] = order + offset + count;
                        ++count;
                }
                offset += count - 1;
        }
}

/*********/
/* TWINS */
/*********/

void TwinReduction::apply() {
        // Vector of hashes of closed neighborhoods
        std::vector<std::pair<size_t, NodeID>> hash_vector;
        hash_vector.reserve(graph_before.number_of_nodes());
        forall_nodes(graph_before, node) {
                hash_vector.push_back({open_neighborhood_hash(graph_before, node), node});
        } endfor

        // Sort by hash
        std::sort(hash_vector.begin(), hash_vector.end());

        // Label nodes to compare neighborhoods
        std::vector<short> labels(graph_before.number_of_nodes(), 0);

        std::vector<std::vector<NodeID>> node_groups;
        node_groups.reserve(graph_before.number_of_nodes());

        // Set to true, if a node has been contracted
        std::vector<bool> contracted(graph_before.number_of_nodes(), false);

        // Store adjusted degrees for each node
        std::vector<NodeWeight> degrees(graph_before.number_of_nodes(), 0);
        forall_nodes(graph_before, node) {
                degrees[node] = compute_reachable_set_size(graph_before, node);
        } endfor

        size_t first = 0, second;
        while (first < hash_vector.size()) {
                auto node_a = hash_vector[first].second;
                second = first + 1;

                node_groups.push_back({node_a});

                forall_out_edges(graph_before, edge, node_a) {
                        labels[graph_before.getEdgeTarget(edge)] = 1;
                } endfor

                while (second < hash_vector.size() && hash_vector[first].first == hash_vector[second].first) {
                        auto node_b = hash_vector[second].second;
                        // ignore contracted nodes
                        if (contracted[node_b]) {
                                second++;
                                continue;
                        }
                        
                        // If neighborhoods are not equal in terms of degrees, we don't need to test
                        if (graph_before.getNodeDegree(node_a) != graph_before.getNodeDegree(node_b) ||
                            degrees[node_a] != degrees[node_b]) {
                                second++;
                                continue;
                        }

                        EdgeWeight count = 0;
                        forall_out_edges(graph_before, edge_b, node_b) {
                                count += labels[graph_before.getEdgeTarget(edge_b)];
                        } endfor
                        if (count == graph_before.getNodeDegree(node_a) && count == graph_before.getNodeDegree(node_b)) {
                                node_groups.back().push_back(node_b);
                                contracted[node_b] = true;
                        }
                        second++;
                }
                do {
                        first++;
                } while(first < hash_vector.size() && contracted[hash_vector[first].second]);
                forall_out_edges(graph_before, edge, node_a) {
                        labels[graph_before.getEdgeTarget(edge)] = 0;
                } endfor
        }

        contract_nodes(graph_before, graph_after, node_groups, mapping);

        // Match reachable set size in reduced and original graph
        for (NodeID new_node_id = 0; new_node_id < node_groups.size(); ++new_node_id) {
                auto reach = degrees[node_groups[new_node_id][0]];
                auto new_reach_without_offset = graph_after.getNodeWeight(new_node_id) - 1;
                forall_out_edges(graph_after, edge, new_node_id) {
                        new_reach_without_offset += graph_after.getNodeWeight(graph_after.getEdgeTarget(edge));
                } endfor
                auto offset = new_reach_without_offset - reach;
                graph_after.set_contraction_offset(new_node_id, offset);
        }
}

void TwinReduction::map(std::vector<NodeID> &reduced_label, std::vector<NodeID> &new_label) const {
        new_label.resize(graph_before.number_of_nodes());
        std::vector<NodeID> orderings(reduced_label.size(), 0);
        for (size_t node = 0; node < reduced_label.size(); node++) {
                orderings[reduced_label[node]] = node;
        }
        NodeID offset = 0;
        for (size_t order = 0; order < orderings.size(); ++order) {
                auto node = orderings[order];
                NodeID count = 0;
                for (const auto original_node: mapping.at(node)) {
                        new_label[original_node] = order + offset + count;
                        ++count;
                }
                offset += count - 1;
        }
}

/********************/
/* PATH COMPRESSION */
/********************/

// Helper functions for finding degree-2 paths
// Recursively walk along a path of degree-2 nodes, pushing each node into 'nodes' along the way

// only detect paths of nodes with weight 1
void degree_2_walk_weightone(graph_access &graph, NodeID start_node, std::vector<NodeID> &nodes) {
        if (graph.getNodeWeight(start_node) != 1) {
                return;
        }
        forall_out_edges(graph, edge, start_node) {
                auto target = graph.getEdgeTarget(edge);
                bool contained = false;
                for (auto node: nodes) {
                        if (target == node) {
                                contained = true;
                        }
                }
                if (!contained && graph.getNodeDegree(target) == 2 && graph.getNodeWeight(target) == 1) {
                        nodes.push_back(target);
                        degree_2_walk_weightone(graph, target, nodes);
                }
        } endfor
}

// detect paths of nodes with any weight
void degree_2_walk_anyweight(graph_access &graph, NodeID start_node, std::vector<NodeID> &nodes) {
        forall_out_edges(graph, edge, start_node) {
                auto target = graph.getEdgeTarget(edge);
                bool contained = false;
                for (auto node: nodes) {
                        if (target == node) {
                                contained = true;
                        }
                }
                if (!contained && graph.getNodeDegree(target) == 2) {
                        nodes.push_back(target);
                        degree_2_walk_anyweight(graph, target, nodes);
                }
        } endfor
}

void PathCompression::apply() {
        std::vector<std::vector<NodeID>> node_groups;
        node_groups.reserve(graph_before.number_of_nodes());
        std::vector<bool> compressed(graph_before.number_of_nodes(), false);

        forall_nodes(graph_before, node) {
                if (graph_before.getNodeDegree(node) != 2) {
                        node_groups.push_back({node});
                } else {
                        if (!compressed[node]) {
                                node_groups.push_back({node});
                                degree_2_walk_weightone(graph_before, node, node_groups.back());
                                for (auto compressed_node: node_groups.back()) {
                                        compressed[compressed_node] = true;
                                }
                        }
                }
        } endfor

        // Contract nodes, setting the weight of new nodes to 1 and their contraction offset to 0
        contract_nodes(graph_before, graph_after, node_groups, mapping);
}

void PathCompression::map(std::vector<NodeID> &reduced_label, std::vector<NodeID> &new_label) const {
        new_label.resize(graph_before.number_of_nodes(), graph_before.number_of_nodes() + 1);
        bucket_pq queue(reduced_label.size());
        for (size_t i = 0; i < reduced_label.size(); ++i) {
                queue.insert(i, -reduced_label[i]);
        }
        NodeID offset = 0;
        while (!queue.empty()) {
                auto order = -queue.maxValue();
                auto node = queue.deleteMax();
                std::vector<NodeID> &path = mapping.at(node);
                if (path.size() == 1) {
                        new_label[path[0]] = order + offset;
                } else {
                        std::array<NodeID, 2> ends;
                        std::array<NodeID, 2> end_neighbors;
                        int endID = 0;
                        int pathID = 0;
                        // find the two end nodes
                        // cycles of degree-2 nodes don't have ends. We take care of that later.
                        for (const auto path_node: path) {
                                forall_out_edges(graph_before, edge, path_node) {
                                        auto target = graph_before.getEdgeTarget(edge);
                                        if (graph_before.getNodeDegree(target) != 2 || graph_before.getNodeWeight(target) != 1) {
                                                end_neighbors[endID] = target;
                                                ends[endID] = pathID;
                                                endID++;
                                        }
                                } endfor
                                pathID++;
                        }
 
                        // If endID > 0, we don't have a cycle of degree-2 nodes
                        // If endID == 0, we have such a cycle and we can just order the nodes arbitrarily
                        if (endID > 0) {
                                // find which non-path neighbor is labeled first,
                                // start ordering from there.
                                // Note that this might still not be optimal, because we don't really know
                                // if any of the path-ends become degree 1. It is possible, that
                                // although end_neighbors[1] is ordered first, ends[0] is simplicial
                                // and we should eliminate from that end.
                                // However, we cannot know that without going through the actual elimination,
                                // so we just try to guess as well as we can.
                                if (new_label.at(end_neighbors.at(0)) < new_label.at(end_neighbors.at(1))) {
                                        endID = 0;
                                } else {
                                        endID = 1;
                                }

                                // reverse the first half of the path, so the nodes are in path order
                                std::reverse(path.begin(), path.begin() + (std::min(ends[0], ends[1]) + 1));

                                // start labeling from the chosen end
                                if (ends[endID] > ends[1 - endID]) {
                                        std::reverse(path.begin(), path.end());
                                }
                        }
                        NodeID count = 0;
                        for (const auto original_node: path) {
                                new_label[original_node] = order + offset + count;
                                ++count;
                        }
                        offset += count - 1;
                }
        }
}

/*************************/
/* DEGREE-2 NODE REMOVAL */
/*************************/

// Find the nodes to replace degree-2 nodes
void find_replacements(graph_access &graph_before, std::vector<std::array<NodeID, 2>> &replacements) {
        replacements.resize(graph_before.number_of_nodes(), {0, 0});
        std::vector<NodeID> path;
        std::vector<bool> visited(graph_before.number_of_nodes(), false);
        forall_nodes(graph_before, node) {
               if(!visited[node] && graph_before.getNodeDegree(node) == 2) {
                        // Detect a path
                        path.clear();
                        path.push_back(node);
                        degree_2_walk_anyweight(graph_before, node, path);
                        // Find the non-path neighbors of the end nodes
                        NodeID ends[2];
                        int i = 0;
                        for (auto path_node: path) {
                                forall_out_edges(graph_before, edge, path_node) {
                                        auto target = graph_before.getEdgeTarget(edge);
                                        if (graph_before.getNodeDegree(target) != 2) {
                                                ends[i] = target;
                                                ++i;
                                        }
                                } endfor
                        }
                        // Set each node on the path to visited and set the found nodes as replacements
                        for (auto path_node: path) {
                                visited[path_node] = true;
                                replacements[path_node][0] = ends[0];
                                replacements[path_node][1] = ends[1];
                        }
                } else if (graph_before.getNodeDegree(node) != 2) {
                        // Make the node its own replacement if it doesn't need to be replaced
                        replacements[node][0] = replacements[node][1] = node;
                }
        } endfor
}

void Degree2Elimination::apply() {
        // Eliminate nodes of degree 2. This is close to what we do in the min-degree algorithm,
        // but no exactly the same:
        //   1. We don't check if we get any nodes of degree 1 in this process. We might lose a bit of quality that way.
        //              If the degree-2 nodes we eliminate are not separators, we don't lose quality at all
        //   2. We don't consider adjusted degree. This is really not a problem though, since
        //      adjusted degree = 2 usually implies degree = 2, except for some cases where degree = 1.
        //      The nodes where degree = 2 and adjusted degree > 2 are contracted indistinguishable nodes.
        //      These nodes are eliminated with fill-in 1, just like regular degree-2 nodes.
        //      Eliminating based on degree instead of adjusted degree also gets rid of more nodes.
        
        graph_after.start_construction(graph_before.number_of_nodes(), graph_before.number_of_edges());
        std::vector<NodeID> reverse_mapping(graph_before.number_of_nodes(), 1000000);
        // Copy nodes of degree != 2, omit nodes of degree 2
        forall_nodes(graph_before, node) {
                if (graph_before.getNodeDegree(node) != 2) {
                        auto new_node_id = graph_after.new_node();
                        graph_after.setNodeWeight(new_node_id, graph_before.getNodeWeight(node));
                        graph_after.set_contraction_offset(new_node_id, graph_before.get_contraction_offset(node));
                        reverse_mapping[node] = new_node_id;
                        mapping.push_back(node);
                } else {
                        label_first.push_back(node);
                }
        } endfor

        std::vector<std::array<NodeID, 2>> replacements;
        find_replacements(graph_before, replacements);

        std::vector<NodeID> update_nodes;
        
        // Copy edges, replacing targets by nodes in the 'replacements' vector if they have degree 2
        // last_visitor stores for a node x the last node y for which an edge (y, x) has been added
        // and the corresponding edge-id
        std::vector<std::pair<NodeID, EdgeID>> last_visitor(graph_before.number_of_nodes(), {-1, -1});
        forall_nodes(graph_before, source) {
                if (graph_before.getNodeDegree(source) == 2) {
                        continue;
                }
                auto new_source_id = reverse_mapping[source];
                bool pushed = false;
                forall_out_edges(graph_before, old_edge, source) {
                        auto old_target = graph_before.getEdgeTarget(old_edge);

                        // If a node has a neighbor that's been eliminated, set its offset to 0.
                        // Its neighborhood is now a clique, due to the elimination process.
                        if (graph_before.getNodeDegree(old_target) == 2) {
                                graph_after.set_contraction_offset(new_source_id, 0);
                                if (!pushed) {
                                        update_nodes.push_back(new_source_id);
                                        pushed = true;
                                }
                        }

                        NodeID new_target_id;
                        // Select the replacement that's not the source
                        if (source != replacements[old_target][0]) {
                                new_target_id = reverse_mapping[replacements[old_target][0]];
                        } else {
                                new_target_id = reverse_mapping[replacements[old_target][1]];
                        }
                        // Add an edge if it's not a self loop and the edge has not been added before
                        if (new_target_id != new_source_id) {
                                if (last_visitor[new_target_id].first != new_source_id) {
                                        auto new_edge_id = graph_after.new_edge(new_source_id, new_target_id);
                                        graph_after.setEdgeWeight(new_edge_id, graph_before.getEdgeWeight(old_edge));
                                        last_visitor[new_target_id] = {new_source_id, new_edge_id};
                                } else {
                                        auto new_edge_id = last_visitor[new_target_id].second;
                                        graph_after.setEdgeWeight(new_edge_id, graph_after.getEdgeWeight(new_edge_id) +
                                                                               graph_before.getEdgeWeight(old_edge));
                                }
                        }
                } endfor
        } endfor

        graph_after.finish_construction();

        // There are cases in which a forward and a backward edge end up with different edge weights:
        // Consider nodes a and b that are connected by two degree-2 nodes x and y.
        // If edges a-x and a-y have weight 1 and edges b-x and b-y have weight 2, the edge a-b will have weight 1
        // and the edge b-a will have weight 2.
        // Identify all these edges and replace their weights by the sum of the weights
        for (const auto node: update_nodes) {
                forall_out_edges(graph_after, forward_edge, node) {
                        auto target = graph_after.getEdgeTarget(forward_edge);
                        forall_out_edges(graph_after, backward_edge, target) {
                                if (graph_after.getEdgeTarget(backward_edge) == node) {
                                        auto forward_weight = graph_after.getEdgeWeight(forward_edge);
                                        auto backward_weight = graph_after.getEdgeWeight(backward_edge);
                                        if (forward_weight != backward_weight) {
                                                graph_after.setEdgeWeight(forward_edge, forward_weight + backward_weight);
                                                graph_after.setEdgeWeight(backward_edge, forward_weight + backward_weight);
                                        }
                                }
                        } endfor
                } endfor
        }
}

void Degree2Elimination::map(std::vector<NodeID> &reduced_label, std::vector<NodeID> &new_label) const {
        new_label.resize(graph_before.number_of_nodes());
        NodeID order = 0;
        for (auto node: label_first) {
                new_label[node] = order;
                order++;
        }
        for (size_t i = 0; i < reduced_label.size(); ++i) {
                new_label[mapping[i]] = reduced_label[i] + order;
        }
}

/************************/
/* TRIANGLE CONTRACTION */
/************************/

// Starting a 'start_node', find neighboring nodes of degree 3 that share at least one neighbor with 'init_node'.
// Node x is ignored if ignore[x] == true
void degree_3_walk(graph_access &graph, NodeID start_node, NodeID init_node, std::vector<NodeID> &path, const std::vector<bool> &ignore) {
        forall_out_edges(graph, edge, start_node) {
                auto target = graph.getEdgeTarget(edge);
                if (ignore[target] || graph.getNodeDegree(target) != 3) {
                        continue;
                }
                bool target_seen = false;
                for (auto n: path) {
                        if (n == target) {
                                target_seen = true;
                                break;
                        }
                }
                if (!target_seen) {
                        // Test if target shares a neighbor with init_node
                        bool shared_neighbor_exists = false;
                        forall_out_edges(graph, edge1, target) {
                                if (!shared_neighbor_exists) {
                                        forall_out_edges(graph, edge2, init_node) {
                                                if (graph.getEdgeTarget(edge1) == graph.getEdgeTarget(edge2)) {
                                                        shared_neighbor_exists = true;
                                                        break;
                                                }
                                        } endfor
                                }
                        } endfor
                        if (shared_neighbor_exists) {
                                // Continue only if there is a shared neighbor 
                                path.push_back(target);
                                degree_3_walk(graph, target, init_node, path, ignore);
                        }
                }
        } endfor
}

void TriangleContraction::apply() {
        std::vector<std::vector<NodeID>> node_groups;
        node_groups.reserve(graph_before.number_of_nodes());
        std::vector<bool> contracted(graph_before.number_of_nodes(), false);
        forall_nodes(graph_before, node) {
                if (contracted[node]) {
                        continue;
                } else if (graph_before.getNodeDegree(node) == 3) {
                        node_groups.push_back({node});
                        degree_3_walk(graph_before, node, node, node_groups.back(), contracted);
                        for (auto n: node_groups.back()) {
                                contracted[n] = true;
                        }
                } else {
                        node_groups.push_back({node});
                }
        } endfor 

        contract_nodes(graph_before, graph_after, node_groups, mapping);
}

void TriangleContraction::map(std::vector<NodeID> &reduced_label, std::vector<NodeID> &new_label) const {
        new_label.resize(graph_before.number_of_nodes());
        std::vector<NodeID> orderings(reduced_label.size(), 0);
        for (size_t node = 0; node < reduced_label.size(); node++) {
                orderings[reduced_label[node]] = node;
        }
        NodeID offset = 0;
        for (size_t order = 0; order < orderings.size(); ++order) {
                auto node = orderings[order];
                NodeID count = 0;
                for (const auto original_node: mapping.at(node)) {
                        new_label[original_node] = order + offset + count;
                        ++count;
                }
                offset += count - 1;
        }
}

/***************************/
/* APPLYING ALL REDUCTIONS */
/***************************/

// Apply the reductions specified in the given configuration to 'in_graph'.
// Store the reductions in 'reductions_stack', in order of application.
// 'recursion_level' is an optional parameter for recording the size of
// reduced graph per recursion level.
void apply_reductions_internal(const PartitionConfig &config,
                               graph_access &in_graph,
                               std::vector<std::unique_ptr<Reduction>> &reduction_stack,
                               int recursion_level = 0) {
        graph_access *graph_1 = &in_graph;              // Point to graph before reduction
        auto num_nodes_before = graph_1->number_of_nodes();
        int round_counter = 0;
        do {
                round_counter++;
                num_nodes_before = graph_1->number_of_nodes();
                for (auto type: config.reduction_order) {
                        timer reduction_timer;
                        // print reduction name
                        switch (type) {
                                case simplicial_nodes:
                                        reduction_stack.emplace_back(new SimplicialNodeReduction(*graph_1, config.max_simplicial_degree));
                                        break;
                                case indistinguishable_nodes:
                                        reduction_stack.emplace_back(new IndistinguishableNodeReduction(*graph_1));
                                        break;
                                case twins:
                                        reduction_stack.emplace_back(new TwinReduction(*graph_1));
                                        break;
                                case path_compression:
                                        reduction_stack.emplace_back(new PathCompression(*graph_1));
                                        break;
                                case degree_2_nodes:
                                        reduction_stack.emplace_back(new Degree2Elimination(*graph_1));
                                        break;
                                case triangle_contraction:
                                        reduction_stack.emplace_back(new TriangleContraction(*graph_1));
                                        break;
                                default:
                                        // Ignore illegal values
                                        continue;
                        }
                        reduction_stack.back()->apply();
                        reduction_stat_counter::get_instance()
                                               .count_reduction(type,
                                                                graph_1->number_of_nodes(),
                                                                reduction_stack.back()->get_reduced_graph().number_of_nodes(),
                                                                recursion_level);

                        graph_1 = &reduction_stack.back()->get_reduced_graph();       
                }
        } while ( (num_nodes_before - reduction_stack.back()->get_reduced_graph().number_of_nodes()) / (double)num_nodes_before >
                                config.convergence_factor);
}

bool apply_reductions(const PartitionConfig &config,
                      graph_access &in_graph,
                      std::vector<std::unique_ptr<Reduction>> &reduction_stack,
                      int recursion_level) {
        if (config.disable_reductions || config.reduction_order.empty()) {
                return false;
        } else {
                apply_reductions_internal(config, in_graph, reduction_stack, recursion_level);
                return true;
        }
}

void map_ordering(const std::vector<std::unique_ptr<Reduction>> &reduction_stack,
                  std::vector<NodeID> &reduced_label,
                  std::vector<NodeID> &final_label) {
        std::vector<NodeID> tmp_label_1, tmp_label_2;
        // Temporary label vectors
        std::vector<NodeID> *label_1, *label_2;
        label_1 = &reduced_label;
        label_2 = &tmp_label_2;
        for (auto reduction = reduction_stack.rbegin(); reduction != reduction_stack.rend(); ++reduction) {
                if (reduction + 1 == reduction_stack.rend()) {
                        label_2 = &final_label;
                }
                (*reduction)->map(*label_1, *label_2);
                if (label_1 == &reduced_label) {
                        label_1 = &tmp_label_1;
                }
                std::swap(label_1, label_2);
        }
}
