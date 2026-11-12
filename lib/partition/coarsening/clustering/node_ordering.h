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
#include "data_structure/graph_access.h"
#include "tools/random_functions.h"
#include "tools/aae_sort.h"

class node_ordering {
public:
        node_ordering();
        virtual ~node_ordering();

        void order_nodes(const PartitionConfig & config, graph_access & G, std::vector< NodeID > & ordered_nodes) {
                forall_nodes(G, node) {
                        ordered_nodes[node] = node;
                } endfor

                switch( config.node_ordering ) {
                        case RANDOM_NODEORDERING:
                                order_nodes_random(config, G, ordered_nodes);
                             break;
                        case DEGREE_NODEORDERING:
                                order_nodes_degree(config, G, ordered_nodes);
                             break;
                 }
        }

        void order_nodes_random(const PartitionConfig & config, graph_access & G, std::vector< NodeID > & ordered_nodes) { 
                random_functions::permutate_vector_fast(ordered_nodes, false);
        }

        void order_nodes_degree(const PartitionConfig & config, graph_access & G, std::vector< NodeID > & ordered_nodes) { 
                // The same std::sort with the same comparisons on (degree, node) entries: std::sort's
                // steps depend only on the comparison results, so the order is the same, but the
                // degrees are read from the entries instead of the node array.
                struct entry {
                        EdgeWeight degree;
                        NodeID     node;
                };
                std::vector<entry> entries(ordered_nodes.size());
                for( size_t i = 0; i < entries.size(); i++) {
                        entries[i].degree = G.getNodeDegree(ordered_nodes[i]);
                        entries[i].node   = ordered_nodes[i];
                }
                // aae_sort: the same arrangement as std::sort, see tools/aae_sort.h
                if(!entries.empty()) aae_sort( entries.data(), entries.data() + entries.size(), 
                           []( const entry & lhs, const entry & rhs) -> bool {
                                return lhs.degree < rhs.degree;
                           });
                for( size_t i = 0; i < entries.size(); i++) ordered_nodes[i] = entries[i].node;
        }

 };


#endif /* end of include guard: NODE_ORDERING_HM1YMLB1 */
