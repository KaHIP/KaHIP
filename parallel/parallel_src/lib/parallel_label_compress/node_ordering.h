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
#include "ppartition_config.h"
#include "data_structure/parallel_graph_access.h"
#include "tools/random_functions.h"

class node_ordering {
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
                std::sort( ordered_nodes.begin(), ordered_nodes.end(), 
                           [&]( const NodeID & lhs, const NodeID & rhs) -> bool {
                                return (G.getNodeDegree(lhs) < G.getNodeDegree(rhs));
                           });
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
