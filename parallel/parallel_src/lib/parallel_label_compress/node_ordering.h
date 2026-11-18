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

// ParHIP's own class, kept apart from the main KaHIP class of the same name
namespace parhip {

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

        // ascending degree, ties in node order (counting sort): nodes of equal degree are visited in
        // memory order, which keeps their own data and edges in cache (std::sort left ties scattered)
        void order_nodes_degree(const PPartitionConfig & config, parallel_graph_access & G, std::vector< NodeID > & ordered_nodes) { 
                NodeID n = ordered_nodes.size();
                EdgeID max_deg = 0;
                for( NodeID node = 0; node < n; node++) max_deg = std::max(max_deg, (EdgeID)G.getNodeDegree(node));
                std::vector< NodeID > start( max_deg + 2, 0 );
                for( NodeID node = 0; node < n; node++) start[G.getNodeDegree(node) + 1]++;
                for( EdgeID d = 1; d < start.size(); d++) start[d] += start[d-1];
                for( NodeID node = 0; node < n; node++) ordered_nodes[start[G.getNodeDegree(node)]++] = node;
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

} // namespace parhip

using parhip::node_ordering;


#endif /* end of include guard: NODE_ORDERING_HM1YMLB1 */
