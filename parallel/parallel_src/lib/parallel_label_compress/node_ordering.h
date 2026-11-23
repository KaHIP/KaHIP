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

        // ascending degree (counting sort); within a degree class, chunks of 128 consecutive nodes in
        // random order: mostly sequential access, but no sweep along the input order
        void order_nodes_degree(const PPartitionConfig & config, parallel_graph_access & G, std::vector< NodeID > & ordered_nodes) { 
                NodeID n = ordered_nodes.size();
                EdgeID max_deg = 0;
                for( NodeID node = 0; node < n; node++) max_deg = std::max(max_deg, (EdgeID)G.getNodeDegree(node));
                std::vector< NodeID > start( max_deg + 2, 0 );
                for( NodeID node = 0; node < n; node++) start[G.getNodeDegree(node) + 1]++;
                for( EdgeID d = 1; d < start.size(); d++) start[d] += start[d-1];
                for( NodeID node = 0; node < n; node++) ordered_nodes[start[G.getNodeDegree(node)]++] = node;

                // start[d] is now the end of class d; shuffle the chunks of each class
                const NodeID chunk = 128;
                std::vector< NodeID > chunks, tmp;
                NodeID begin = 0;
                for( EdgeID d = 0; d <= max_deg; d++) {
                        NodeID end = start[d];
                        NodeID num_chunks = (end - begin + chunk - 1) / chunk;
                        if( num_chunks > 1 ) {
                                chunks.resize(num_chunks);
                                for( NodeID c = 0; c < num_chunks; c++) chunks[c] = c;
                                for( NodeID c = num_chunks - 1; c > 0; c--) std::swap(chunks[c], chunks[random_functions::nextInt<NodeID>(0, c)]);
                                tmp.assign(ordered_nodes.begin() + begin, ordered_nodes.begin() + end);
                                NodeID pos = begin;
                                for( NodeID c = 0; c < num_chunks; c++) {
                                        NodeID from = chunks[c] * chunk, to = std::min(from + chunk, end - begin);
                                        for( NodeID x = from; x < to; x++) ordered_nodes[pos++] = tmp[x];
                                }
                        }
                        begin = end;
                }
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
