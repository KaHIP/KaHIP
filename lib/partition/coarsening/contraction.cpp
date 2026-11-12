/******************************************************************************
 * contraction.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <memory>

#include "contraction.h"
#include "../uncoarsening/refinement/quotient_graph_refinement/complete_boundary.h"
#include "macros_assertions.h"

contraction::contraction() {

}

contraction::~contraction() {

}

// for documentation see technical reports of christian schulz  
void contraction::contract(const PartitionConfig & partition_config, 
                           graph_access & G, 
                           graph_access & coarser, 
                           const Matching & edge_matching,
                           const CoarseMapping & coarse_mapping,
                           const NodeID & no_of_coarse_vertices,
                           const NodePermutationMap & permutation) const {

        if(partition_config.matching_type == CLUSTER_COARSENING) {
                return contract_clustering(partition_config, G, coarser, edge_matching, coarse_mapping, no_of_coarse_vertices, permutation);
        }

        if(partition_config.combine) {
                coarser.resizeSecondPartitionIndex(no_of_coarse_vertices);
        }

        // targets in the coarse graph per edge; every entry is written below, so the buffer is not
        // value-initialized
        const EdgeID m = G.number_of_edges();
        std::unique_ptr<NodeID[]> new_edge_targets_buf(new NodeID[m]);
        NodeID * new_edge_targets = new_edge_targets_buf.get();
        forall_edges(G, e) {
                new_edge_targets[e] = coarse_mapping[G.getEdgeTarget(e)];
        } endfor

        std::vector<EdgeID> edge_positions(no_of_coarse_vertices, UNDEFINED_EDGE);

        // The coarse edges are collected first (coarse node c owns c_edges[c_first[c] .. c_first[c+1])),
        // so the coarse graph is built with its real edge count instead of allocating and zeroing m
        // edges and m ratings. Same nodes, edges, order and weights as before.
        struct c_edge {
                NodeID     target;
                EdgeWeight weight;
        };
        std::vector<c_edge> c_edges;
        c_edges.reserve(m);
        std::vector<EdgeID>     c_first(no_of_coarse_vertices + 1, 0);
        std::vector<NodeWeight> c_weight(no_of_coarse_vertices);

        NodeID cur_no_vertices = 0;
        auto visit = [&](NodeID coarseNode, EdgeID e) {
                const NodeID t = new_edge_targets[e];
                if(t == coarseNode) return; //this is the matched edge ... return
                EdgeID & edge_pos = edge_positions[t];
                if( edge_pos == UNDEFINED_EDGE ) {
                        edge_pos = c_edges.size();
                        c_edges.push_back(c_edge{t, G.getEdgeWeight(e)});
                } else {
                        c_edges[edge_pos].weight += G.getEdgeWeight(e);
                }
        };

        forall_nodes(G, n) {
                NodeID node = permutation[n];
                //we look only at the coarser nodes
                if(coarse_mapping[node] != cur_no_vertices) 
                        continue;
                
                NodeID coarseNode = cur_no_vertices;
                c_weight[coarseNode] = G.getNodeWeight(node);

                if(partition_config.combine) {
                        coarser.setSecondPartitionIndex(coarseNode, G.getSecondPartitionIndex(node));
                }

                // do something with all outgoing edges (in auxillary graph)
                forall_out_edges(G, e, node) {
                        visit(coarseNode, e);
                } endfor

                //this node was really matched
                NodeID matched_neighbor = edge_matching[node];
                if(node != matched_neighbor) {
                        //update weight of coarser node
                        c_weight[coarseNode] = G.getNodeWeight(node) + G.getNodeWeight(matched_neighbor);

                        forall_out_edges(G, e, matched_neighbor) {
                                visit(coarseNode, e);
                        } endfor
                }
                for( EdgeID i = c_first[coarseNode]; i < c_edges.size(); i++) {
                       edge_positions[c_edges[i].target] = UNDEFINED_EDGE;
                }
                
                cur_no_vertices++;
                c_first[cur_no_vertices] = c_edges.size();
        } endfor

        ASSERT_RANGE_EQ(edge_positions, 0, edge_positions.size(), UNDEFINED_EDGE); 
        ASSERT_EQ(no_of_coarse_vertices, cur_no_vertices);

        coarser.start_construction(no_of_coarse_vertices, c_edges.size());
        for( NodeID c = 0; c < no_of_coarse_vertices; c++) {
                NodeID coarseNode = coarser.new_node();
                coarser.setNodeWeight(coarseNode, c_weight[c]);
                for( EdgeID i = c_first[c]; i < c_first[c+1]; i++) {
                        EdgeID coarse_edge = coarser.new_edge(coarseNode, c_edges[i].target);
                        coarser.setEdgeWeight(coarse_edge, c_edges[i].weight);
                }
        }
        
        //this also resizes the edge fields ... 
        coarser.finish_construction();
}

void contraction::contract_clustering(const PartitionConfig & partition_config, 
                              graph_access & G, 
                              graph_access & coarser, 
                              const Matching & edge_matching,
                              const CoarseMapping & coarse_mapping,
                              const NodeID & no_of_coarse_vertices,
                              const NodePermutationMap & permutation) const {

        if(partition_config.combine) {
                coarser.resizeSecondPartitionIndex(no_of_coarse_vertices);
        }

        //save partition map -- important if the graph is allready partitioned
        std::vector< int > partition_map(G.number_of_nodes());
        int k = G.get_partition_count();
        forall_nodes(G, node) {
                partition_map[node] = G.getPartitionIndex(node);
                G.setPartitionIndex(node, coarse_mapping[node]);
        } endfor

        G.set_partition_count(no_of_coarse_vertices);

        complete_boundary bnd(&G);
        //bnd.build();
        bnd.fastComputeQuotientGraph(coarser, no_of_coarse_vertices);
        //bnd.getUnderlyingQuotientGraph(coarser);
        //

        G.set_partition_count(k);
        forall_nodes(G, node) {
                G.setPartitionIndex(node, partition_map[node]);
                coarser.setPartitionIndex(coarse_mapping[node], G.getPartitionIndex(node));

                if(partition_config.combine) {
                        coarser.setSecondPartitionIndex(coarse_mapping[node], G.getSecondPartitionIndex(node));
                }

        } endfor

}


// for documentation see technical reports of christian schulz  
void contraction::contract_partitioned(const PartitionConfig & partition_config, 
                                       graph_access & G, 
                                       graph_access & coarser, 
                                       const Matching & edge_matching,
                                       const CoarseMapping & coarse_mapping,
                                       const NodeID & no_of_coarse_vertices,
                                       const NodePermutationMap & permutation) const {
        
        if(partition_config.matching_type == CLUSTER_COARSENING) {
                return contract_clustering(partition_config, G, coarser, edge_matching, coarse_mapping, no_of_coarse_vertices, permutation);
        }


        std::vector<NodeID> new_edge_targets(G.number_of_edges());
        forall_edges(G, e) {
                new_edge_targets[e] = coarse_mapping[G.getEdgeTarget(e)];
        } endfor

        std::vector<EdgeID> edge_positions(no_of_coarse_vertices, UNDEFINED_EDGE);

        //we dont know the number of edges jet, so we use the old number for 
        //construction of the coarser graph and then resize the field according
        //to the number of edges we really got
        coarser.set_partition_count(G.get_partition_count());
        coarser.start_construction(no_of_coarse_vertices, G.number_of_edges());

        if(partition_config.combine) {
                coarser.resizeSecondPartitionIndex(no_of_coarse_vertices);
        }

        NodeID cur_no_vertices = 0;

        PRINT(std::cout <<  "contracting a partitioned graph"  << std::endl;)
        forall_nodes(G, n) {
                NodeID node = permutation[n];
                //we look only at the coarser nodes
                if(coarse_mapping[node] != cur_no_vertices) 
                        continue;
                
                NodeID coarseNode = coarser.new_node();
                coarser.setNodeWeight(coarseNode, G.getNodeWeight(node));
                coarser.setPartitionIndex(coarseNode, G.getPartitionIndex(node));

                if(partition_config.combine) {
                        coarser.setSecondPartitionIndex(coarseNode, G.getSecondPartitionIndex(node));
                }
                // do something with all outgoing edges (in auxillary graph)
                forall_out_edges(G, e, node) {
                                visit_edge(G, coarser, edge_positions, coarseNode, e, new_edge_targets);                        
                } endfor

                //this node was really matched
                NodeID matched_neighbor = edge_matching[node];
                if(node != matched_neighbor) {
                        //update weight of coarser node
                        NodeWeight new_coarse_weight = G.getNodeWeight(node) + G.getNodeWeight(matched_neighbor);
                        coarser.setNodeWeight(coarseNode, new_coarse_weight);

                        forall_out_edges(G, e, matched_neighbor) {
                                visit_edge(G, coarser, edge_positions, coarseNode, e, new_edge_targets);
                        } endfor
                }
                forall_out_edges(coarser, e, coarseNode) {
                       edge_positions[coarser.getEdgeTarget(e)] = UNDEFINED_EDGE;
                } endfor
                
                cur_no_vertices++;
        } endfor

        ASSERT_RANGE_EQ(edge_positions, 0, edge_positions.size(), UNDEFINED_EDGE); 
        ASSERT_EQ(no_of_coarse_vertices, cur_no_vertices);
        
        //this also resizes the edge fields ... 
        coarser.finish_construction();
}



