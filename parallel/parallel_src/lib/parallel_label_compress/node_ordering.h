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

        // std::sort of the nodes by degree. The degree is packed into the upper 32 bits of a key
        // with the node in the lower bits and only the upper bits are compared, so std::sort sees
        // the same comparison results as with the original comparator on the node ids and
        // produces the same permutation, without loading two degrees per comparison.
        void sort_by_degree(parallel_graph_access & G, std::vector< NodeID > & ordered_nodes) {
                std::vector< unsigned long long > keys( ordered_nodes.size() );
                for( size_t i = 0; i < ordered_nodes.size(); i++) {
                        keys[i] = ((unsigned long long) G.getNodeDegree(ordered_nodes[i]) << 32) | ordered_nodes[i];
                }
                std::sort( keys.begin(), keys.end(), 
                           []( const unsigned long long & lhs, const unsigned long long & rhs) -> bool {
                                return (lhs >> 32) < (rhs >> 32);
                           });
                for( size_t i = 0; i < ordered_nodes.size(); i++) {
                        ordered_nodes[i] = (NodeID)(keys[i] & 0xFFFFFFFFULL);
                }
        }

        void order_nodes_random(const PPartitionConfig & config, parallel_graph_access & G, std::vector< NodeID > & ordered_nodes) { 
                random_functions::permutate_vector_fast(ordered_nodes, false);
        }

        void order_nodes_degree(const PPartitionConfig & config, parallel_graph_access & G, std::vector< NodeID > & ordered_nodes) { 
                sort_by_degree(G, ordered_nodes);
        }

        void order_leastghostnodes_nodes_degree(const PPartitionConfig & config, parallel_graph_access & G, std::vector< NodeID > & ordered_nodes) { 
                sort_by_degree(G, ordered_nodes);
        }

	void order_nodes_degree_leastghostnodes(const PPartitionConfig & config, parallel_graph_access & G, std::vector< NodeID > & ordered_nodes) { 
                sort_by_degree(G, ordered_nodes);
        }
};

} // namespace parhip

using parhip::node_ordering;


#endif /* end of include guard: NODE_ORDERING_HM1YMLB1 */
