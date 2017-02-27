/******************************************************************************
 * contraction.cpp
 *
 * Source of KaHIP -- Karlsruhe High Quality Graph Partitioning 
 ******************************************************************************
 * Copyright (C) 2017 Christian Schulz 
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

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

        std::vector<NodeID> new_edge_targets(G.number_of_edges());
        forall_edges(G, e) {
                new_edge_targets[e] = coarse_mapping[G.getEdgeTarget(e)];
        } endfor

        std::vector<EdgeID> edge_positions(no_of_coarse_vertices, UNDEFINED_EDGE);

        //we dont know the number of edges jet, so we use the old number for 
        //construction of the coarser graph and then resize the field according
        //to the number of edges we really got
        coarser.start_construction(no_of_coarse_vertices, G.number_of_edges());

        NodeID cur_no_vertices = 0;

        forall_nodes(G, n) {
                NodeID node = permutation[n];
                //we look only at the coarser nodes
                if(coarse_mapping[node] != cur_no_vertices) 
                        continue;
                
                NodeID coarseNode = coarser.new_node();
                coarser.setNodeWeight(coarseNode, G.getNodeWeight(node));

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
        bnd.build();
        bnd.getUnderlyingQuotientGraph(coarser);

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



