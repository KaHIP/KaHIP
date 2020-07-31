/*
 * Author: Wolfgang Ost
 */

#ifndef MIN_DEGREE_ORDERING
#define MIN_DEGREE_ORDERING

#include <unordered_set>
#include <vector>

#include "data_structure/graph_access.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "definitions.h"

// A clique, described by its set of nodes
class clique {

public:
        // create a clique from a list of nodes
        inline clique(std::unordered_set<NodeID> n = {}) : nodes(n) {};

        // check if a node is part of this clique
        inline bool contains(NodeID node) const { 
                return nodes.find(node) != nodes.end();
        }

        template<class InputIt>
        inline bool contains_all_of(InputIt first, InputIt last) {
                for (auto it = first; it != last; ++it) {
                        if (!contains(*it)) {
                                return false;
                        }
                }
                return true;
        }

        // merge the nodes of other into this clique
        inline void merge(clique other) {
                nodes.insert(other.nodes.begin(), other.nodes.end());
        }

        // remove node from this clique
        inline std::unordered_set<NodeID>::const_iterator remove(NodeID node) {
                auto it = nodes.find(node);
                if (it != nodes.end()) {
                        return nodes.erase(it);
                }
                return nodes.end();
        }

        // add node to this clique
        inline void insert(NodeID node) {
                nodes.insert(node);
        }

        inline bool empty() const {
                return nodes.empty();
        }

        // return the number of nodes in the clique
        inline std::unordered_set<NodeID>::size_type size() const {
                return nodes.size();
        }

        inline bool operator==(const clique &other) const {
                return nodes == other.nodes;
        }

        // iterators over the nodes in the clique
        inline const std::unordered_set<NodeID>::const_iterator begin() const {
                return nodes.cbegin();
        }

        inline const std::unordered_set<NodeID>::const_iterator end() const {
                return nodes.cend();
        }

private:
        std::unordered_set<NodeID> nodes;

};


// Use the min degree algorithm to compute a fill-reducing node ordering of G
// Expectes 'labels' to have the correct size (= number of nodes).
// Input:  G, the graph to be ordered
// Output: labels, the order assigned to the nodes
void min_degree_ordering(graph_access &G, std::vector<NodeID> &labels);

// Reimplementation of the minimum degree algorithm
// to get a min-degree ordering for some graph G, call 'MinDegree(&G).perform_ordering(labels);'
//
// This implements:
//  - element absorption
//  - a simplified version of simplicial node reduction
//  - indistinguishable nodes, following George, Liu 1980, "A Fast Implementation of the Minimum Degree Algorithm Using Quotient Graphs"
class MinDegree {

public:
        MinDegree(graph_access * const graph, const std::vector<NodeID> &halo_nodes = {});

        // Order the graph. The result is stored in 'labels'.
        // This vector needs to be set to the number of nodes in the graph, including halo nodes.
        void perform_ordering(std::vector<NodeID> &labels);

private:
        graph_access * const graph;

        // Nodes that haven't been eliminated
        // This is mainly used to eliminate nodes that are only members of one clique
        std::unordered_set<NodeID> current_nodes;

        // This "clique" tracks all neighbors of eliminated nodes that are still in the graph
        // It's not actually a clique, but the clique-class is useful here
        clique elimination_neighbors;

        // Priority queue to store node degrees
        maxNodeHeap degree_pq;

        // Cliques used to represent the graph. Cliques never get deleted,
        // but cliques that aren't used are referenced in the memberships vector
        std::vector<clique> cliques;
        // For each node, memberships stores the indices of cliques the node belongs to
        std::vector<std::vector<size_t>> memberships;

        // For each node, store a link to an indistinguishable node it represents, or itself.
        std::vector<NodeID> node_links;
        // Sum of the weights of indistinguishable nodes, per representative
        std::vector<NodeWeight> link_weights;
        std::vector<NodeID> new_indistinguishable_nodes;

        std::vector<bool> is_halo_node;

        // Compute the initial adjusted degrees for the given graph
        void initialize_degrees();
        // Create the initial clique cover
        void initialize_cliques();

        // Eliminate a node
        // This merges all cliques containing node
        // Returns the index of the merged clique
        size_t eliminate_node(NodeID node);

        // update node degrees of neighbors of the eliminated node.
        // All these neighbors are in the merged clique.
        // This also detects indistinguishable nodes and stores them in new_indistinguishable_nodes.
        void update_node_degrees(NodeID eliminatedNode, size_t merged_clique_idx);

        // Taking the nodes found to be indistinguishable ('new_indistinguishable_nodes'),
        // update the vectors 'node_links' and 'link_weights',
        // and remove the nodes from the relevant sets.
        void update_indistinguishable_nodes();

        // Label the given node and all its indistinguishable nodes
        void label_node(NodeID node, std::vector<NodeID> &labels, NodeID &order);

};

#endif
