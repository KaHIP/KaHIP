/******************************************************************************
 * random_matching.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "macros_assertions.h"
#include "random_functions.h"
#include "random_matching.h"

random_matching::random_matching() {

}

random_matching::~random_matching() {

}

void random_matching::match(const PartitionConfig & partition_config, 
                            graph_access & G, 
                            Matching & edge_matching,
                            CoarseMapping & coarse_mapping,
                            NodeID & no_of_coarse_vertices,
                            NodePermutationMap & permutation) {


        permutation.resize(G.number_of_nodes());
        edge_matching.resize(G.number_of_nodes());
        coarse_mapping.resize(G.number_of_nodes());

        no_of_coarse_vertices = 0;

        if(!(partition_config.matching_type == MATCHING_RANDOM_GPA)) { 
                random_functions::permutate_entries(partition_config, permutation, true);
        } else {
                for( unsigned int i = 0; i < permutation.size(); i++) {
                        permutation[i] = i;
                }                   
        }

        forall_nodes(G, n) {
                edge_matching[n] = n;
        } endfor

        if(partition_config.graph_allready_partitioned) { //in this case edges between partitions arent matched 
                forall_nodes(G, n) {
                        NodeID curNode = permutation[n];
                        NodeWeight curNodeWeight = G.getNodeWeight(curNode);

                        if(edge_matching[curNode] == curNode) {
                                //match with a random neighbor
                                int matchingPartner = curNode;
                                forall_out_edges(G, e, curNode) {
                                        NodeID target = G.getEdgeTarget(e);
                                        NodeWeight coarser_weight = G.getNodeWeight(target) + curNodeWeight;

                                        if(edge_matching[target] == target 
                                        && coarser_weight <= partition_config.max_vertex_weight) {
                                                if(G.getPartitionIndex(curNode) != G.getPartitionIndex(target))
                                                        continue;

                                                if(partition_config.combine) {
                                                        if(G.getSecondPartitionIndex(curNode) != G.getSecondPartitionIndex(target))
                                                                continue;
                                                }

                                                matchingPartner = target;
                                                ASSERT_NEQ(curNode, target);
                                                break;
                                        }
                                } endfor

                                coarse_mapping[matchingPartner] = no_of_coarse_vertices;
                                coarse_mapping[curNode]         = no_of_coarse_vertices;

                                edge_matching[matchingPartner] = curNode;
                                edge_matching[curNode]         = matchingPartner;

                                no_of_coarse_vertices++;
                        } 
                } endfor
        } else {
                //copy n paste from the first if clause but this time all edges are matchable
                forall_nodes(G, n) {
                        NodeID curNode = permutation[n];
                        NodeWeight curNodeWeight = G.getNodeWeight(curNode);

                        if(edge_matching[curNode] == curNode) {
                                //match with a random neighbor
                                int matchingPartner = curNode;
                                forall_out_edges(G, e, curNode) {
                                        NodeID target             = G.getEdgeTarget(e);
                                        NodeWeight coarser_weight = G.getNodeWeight(target) + curNodeWeight;

                                        if(edge_matching[target] == target 
                                        && coarser_weight <= partition_config.max_vertex_weight) {
                                                matchingPartner = target;
                                                ASSERT_NEQ(curNode, target);
                                                break;
                                        }
                                } endfor

                                coarse_mapping[matchingPartner] = no_of_coarse_vertices;
                                coarse_mapping[curNode]         = no_of_coarse_vertices;

                                edge_matching[matchingPartner] = curNode;
                                edge_matching[curNode]         = matchingPartner;

                                no_of_coarse_vertices++;
                        } 
                } endfor

        }
        PRINT(std::cout << "log>" << "no of coarse nodes: " << no_of_coarse_vertices << std::endl;)
}
