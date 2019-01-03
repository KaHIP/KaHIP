/******************************************************************************
 * graph_extractor.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <unordered_map>
#include "graph_extractor.h"


graph_extractor::graph_extractor() {

}

graph_extractor::~graph_extractor() {

}

void graph_extractor::extract_block(graph_access & G, 
                                    graph_access & extracted_block, 
                                    PartitionID block, 
                                    std::vector<NodeID> & mapping) {

        // build reverse mapping
        std::vector<NodeID> reverse_mapping;
        NodeID nodes = 0;
        NodeID dummy_node = G.number_of_nodes() + 1;
        forall_nodes(G, node) {
                if(G.getPartitionIndex(node) == block) {
                        reverse_mapping.push_back(nodes++);
                } else {
                        reverse_mapping.push_back(dummy_node);
                }
        } endfor

        extracted_block.start_construction(nodes, G.number_of_edges());

        forall_nodes(G, node) {
                if(G.getPartitionIndex(node) == block) {
                        NodeID new_node = extracted_block.new_node();
                        mapping.push_back(node);
                        extracted_block.setNodeWeight( new_node, G.getNodeWeight(node));

                        forall_out_edges(G, e, node) {
                                NodeID target = G.getEdgeTarget(e);
                                if( G.getPartitionIndex( target ) == block ) {
                                        EdgeID new_edge = extracted_block.new_edge(new_node, reverse_mapping[target]);
                                        extracted_block.setEdgeWeight(new_edge, G.getEdgeWeight(e));
                                }
                        } endfor
                }
        } endfor

        extracted_block.finish_construction();
}


void graph_extractor::extract_two_blocks(graph_access & G, 
                                         graph_access & extracted_block_lhs, 
                                         graph_access & extracted_block_rhs, 
                                         std::vector<NodeID> & mapping_lhs,
                                         std::vector<NodeID> & mapping_rhs,
                                         NodeWeight & partition_weight_lhs,
                                         NodeWeight & partition_weight_rhs) {

        PartitionID lhs = 0;
        PartitionID rhs = 1;

        // build reverse mapping
        std::vector<NodeID> reverse_mapping_lhs;
        std::vector<NodeID> reverse_mapping_rhs;
        NodeID nodes_lhs     = 0;
        NodeID nodes_rhs     = 0;
        partition_weight_lhs = 0;
        partition_weight_rhs = 0;
        NodeID dummy_node    = G.number_of_nodes() + 1;

        forall_nodes(G, node) {
                if(G.getPartitionIndex(node) == lhs) {
                        reverse_mapping_lhs.push_back(nodes_lhs++);
                        reverse_mapping_rhs.push_back(dummy_node);
                        partition_weight_lhs += G.getNodeWeight(node);
                } else {
                        reverse_mapping_rhs.push_back(nodes_rhs++);
                        reverse_mapping_lhs.push_back(dummy_node);
                        partition_weight_rhs += G.getNodeWeight(node);
                }
        } endfor

        extracted_block_lhs.start_construction(nodes_lhs, G.number_of_edges());
        extracted_block_rhs.start_construction(nodes_rhs, G.number_of_edges());

        forall_nodes(G, node) {
                if(G.getPartitionIndex(node) == lhs) {
                        NodeID new_node = extracted_block_lhs.new_node();
                        mapping_lhs.push_back(node);
                        extracted_block_lhs.setNodeWeight(new_node, G.getNodeWeight(node));

                        forall_out_edges(G, e, node) {
                                NodeID target = G.getEdgeTarget(e);
                                if( G.getPartitionIndex( target ) == lhs) {
                                        EdgeID new_edge = extracted_block_lhs.new_edge(new_node, reverse_mapping_lhs[target]);
                                        extracted_block_lhs.setEdgeWeight( new_edge, G.getEdgeWeight(e));
                                }
                        } endfor

                } else {
                        NodeID new_node = extracted_block_rhs.new_node();
                        mapping_rhs.push_back(node);
                        extracted_block_rhs.setNodeWeight(new_node, G.getNodeWeight(node));

                        forall_out_edges(G, e, node) {
                                NodeID target = G.getEdgeTarget(e);
                                if( G.getPartitionIndex( target ) == rhs) {
                                        EdgeID new_edge = extracted_block_rhs.new_edge(new_node, reverse_mapping_rhs[target]);
                                        extracted_block_rhs.setEdgeWeight( new_edge, G.getEdgeWeight(e));
                                }
                        } endfor
                }
        } endfor

        extracted_block_lhs.finish_construction();
        extracted_block_rhs.finish_construction();
}

// Method takes a number of nodes and extracts the underlying subgraph from G
// it also assignes block informations
void graph_extractor::extract_two_blocks_connected(graph_access & G, 
                                                   std::vector<NodeID> lhs_nodes,
                                                   std::vector<NodeID> rhs_nodes,
                                                   PartitionID lhs, 
                                                   PartitionID rhs,
                                                   graph_access & pair,
                                                   std::vector<NodeID> & mapping) {
        //// build reverse mapping
        std::unordered_map<NodeID,NodeID> reverse_mapping;
        NodeID nodes = 0;
        EdgeID edges = 0; // upper bound for number of edges

        for( unsigned i = 0; i < lhs_nodes.size(); i++) {
                NodeID node           = lhs_nodes[i];
                reverse_mapping[node] = nodes;
                edges += G.getNodeDegree(lhs_nodes[i]);
                nodes++;
        }
        for( unsigned i = 0; i < rhs_nodes.size(); i++) {
                NodeID node           = rhs_nodes[i];
                reverse_mapping[node] = nodes;
                edges += G.getNodeDegree(rhs_nodes[i]);
                nodes++;
        }

        pair.start_construction(nodes, edges);

        for( unsigned i = 0; i < lhs_nodes.size(); i++) {
                NodeID node     = lhs_nodes[i];
                NodeID new_node = pair.new_node();
                mapping.push_back(node);

                pair.setNodeWeight(new_node, G.getNodeWeight(node));
                pair.setPartitionIndex(new_node, 0);

                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if( G.getPartitionIndex( target ) == lhs || G.getPartitionIndex( target ) == rhs ) {
                                EdgeID new_edge = pair.new_edge(new_node, reverse_mapping[target]);
                                pair.setEdgeWeight(new_edge, G.getEdgeWeight(e));
                        }
                } endfor

        }

        for( unsigned i = 0; i < rhs_nodes.size(); i++) {
                NodeID node     = rhs_nodes[i];
                NodeID new_node = pair.new_node();
                mapping.push_back(node);

                pair.setNodeWeight(new_node, G.getNodeWeight(node));
                pair.setPartitionIndex(new_node, 1);

                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if( G.getPartitionIndex( target ) == lhs || G.getPartitionIndex( target ) == rhs ) {
                                EdgeID new_edge = pair.new_edge(new_node, reverse_mapping[target]);
                                pair.setEdgeWeight(new_edge, G.getEdgeWeight(e));
                        }
                } endfor

        }

        pair.finish_construction();
}


