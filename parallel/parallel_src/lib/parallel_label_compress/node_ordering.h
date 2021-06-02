/******************************************************************************
 * node_ordering.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef NODE_ORDERING_HM1YMLB1
#define NODE_ORDERING_HM1YMLB1

#include <algorithm>

#include "definitions.h"
#include "partition_config.h"
#include "data_structure/parallel_graph_access.h"
#include "tools/random_functions.h"

class node_ordering {
private:
        void order_nodes_degree_with_counting(const PPartitionConfig & config, parallel_graph_access & G, std::vector< NodeID > & ordered_nodes) {
                if (ordered_nodes.empty()) return;

                EdgeID offset = std::numeric_limits< EdgeID >::max();
                size_t count = ordered_nodes.size();

                for (NodeID node : ordered_nodes) {
                        EdgeID degree = G.getNodeDegree(node);
                        if (degree < offset) {
                                offset = degree;
                        }
                }
                
                std::vector< size_t > counting_sort_count;
                std::vector< NodeID > node_buffer;

                // last bucket is for nodes with large degree, and will be sorted at the end with std::sort
                counting_sort_count.assign(count + 1, 0);
                for (NodeID node : ordered_nodes) {
                        EdgeID degree = G.getNodeDegree(node);
                        size_t scaled_degree = std::min<size_t>(degree - offset, count);
                        ++counting_sort_count[scaled_degree];
                }
                for (int i = 1; i <= count; i++) {
                        counting_sort_count[i] += counting_sort_count[i - 1];
                }

                node_buffer.resize(ordered_nodes.size());
                for (NodeID node : ordered_nodes) {
                        EdgeID degree = G.getNodeDegree(node);
                        size_t scaled_degree = std::min<size_t>(degree - offset, count);
                        size_t pos = --counting_sort_count[scaled_degree];
                        node_buffer[pos] = node;
                }
                ordered_nodes.assign(node_buffer.begin(), node_buffer.end());
                     
                std::sort( ordered_nodes.begin() + counting_sort_count.back(), ordered_nodes.end(), 
                           [&]( const NodeID & lhs, const NodeID & rhs) -> bool {
                                return (G.getNodeDegree(lhs) < G.getNodeDegree(rhs));
                           });
        }
public:
        node_ordering();
        virtual ~node_ordering();

        void order_nodes(const PPartitionConfig & config, parallel_graph_access & G, std::vector< NodeID > & ordered_nodes) {
                forall_local_nodes(G, node) {
                        ordered_nodes[node] = node;
                } endfor

                switch( config.node_ordering ) {
                        case RANDOM_NODEORDERING:
                                order_nodes_random(config, G, ordered_nodes);
                             break;
                        case DEGREE_NODEORDERING:
                                order_nodes_degree(config, G, ordered_nodes);
                             break;
                        case LEASTGHOSTNODESFIRST_DEGREE_NODEODERING:
                                order_leastghostnodes_nodes_degree(config, G, ordered_nodes);
                             break;
                        case DEGREE_LEASTGHOSTNODESFIRST_NODEODERING:
                                order_nodes_degree_leastghostnodes(config, G, ordered_nodes);
                             break;
                 }
        }

        void order_nodes_random(const PPartitionConfig & config, parallel_graph_access & G, std::vector< NodeID > & ordered_nodes) { 
                random_functions::permutate_vector_fast(ordered_nodes, false);
        }

        void order_nodes_degree(const PPartitionConfig & config, parallel_graph_access & G, std::vector< NodeID > & ordered_nodes) { 
                order_nodes_degree_with_counting(config, G, ordered_nodes);
        }

        void order_leastghostnodes_nodes_degree(const PPartitionConfig & config, parallel_graph_access & G, std::vector< NodeID > & ordered_nodes) { 
                std::sort( ordered_nodes.begin(), ordered_nodes.end(), 
                           [&]( const NodeID & lhs, const NodeID & rhs) -> bool {
                                return (G.getNodeDegree(lhs) < G.getNodeDegree(rhs));
                           });
        }

	void order_nodes_degree_leastghostnodes(const PPartitionConfig & config, parallel_graph_access & G, std::vector< NodeID > & ordered_nodes) { 
                std::sort( ordered_nodes.begin(), ordered_nodes.end(), 
                           [&]( const NodeID & lhs, const NodeID & rhs) -> bool {
				return (G.getNodeDegree(lhs) < G.getNodeDegree(rhs));
                           });
        }
};


#endif /* end of include guard: NODE_ORDERING_HM1YMLB1 */
