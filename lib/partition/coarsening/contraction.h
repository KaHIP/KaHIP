/******************************************************************************
 * contraction.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef CONTRACTION_VIXZ9K0F
#define CONTRACTION_VIXZ9K0F


#include "data_structure/graph_access.h"
#include "matching/matching.h"
#include "partition_config.h"

typedef NodeID Regions;

class contraction {
        public:
                contraction();
                virtual ~contraction();

                void contract(const PartitionConfig & partition_config, 
                              graph_access & finer, 
                              graph_access & coarser, 
                              const Matching & edge_matching,
                              const CoarseMapping & coarse_mapping,
                              const NodeID & no_of_coarse_vertices,
                              const NodePermutationMap & permutation) const;

                void contract_clustering(const PartitionConfig & partition_config, 
                              graph_access & finer, 
                              graph_access & coarser, 
                              const Matching & edge_matching,
                              const CoarseMapping & coarse_mapping,
                              const NodeID & no_of_coarse_vertices,
                              const NodePermutationMap & permutation) const;


                 void contract_partitioned(const PartitionConfig & partition_config, 
                                           graph_access & G, 
                                           graph_access & coarser, 
                                           const Matching & edge_matching,
                                           const CoarseMapping & coarse_mapping,
                                           const NodeID & no_of_coarse_vertices,
                                           const NodePermutationMap & permutation) const; 

        private:
                // visits an edge in G (and auxillary graph) and updates/creates and edge in coarser graph 
                void visit_edge(graph_access & G, 
                                graph_access & coarser,
                                std::vector<EdgeID> & edge_positions,
                                const NodeID coarseNode,
                                const EdgeID e,
                                const std::vector<NodeID> & new_edge_targets) const;


};

inline void contraction::visit_edge(graph_access & G, 
                graph_access & coarser,
                std::vector<EdgeID> & edge_positions,
                const NodeID coarseNode,
                const EdgeID e,
                const std::vector<NodeID> & new_edge_targets) const {

        EdgeID new_coarse_edge_target = new_edge_targets[e];
        if(new_coarse_edge_target == coarseNode) return; //this is the matched edge ... return

        EdgeID edge_pos = edge_positions[new_coarse_edge_target];
        if( edge_pos == UNDEFINED_EDGE ) {
                //we havent seen this target node before so we need to create an edge
                EdgeID coarse_edge = coarser.new_edge(coarseNode, new_coarse_edge_target);
                coarser.setEdgeWeight(coarse_edge, G.getEdgeWeight(e));
                edge_positions[new_coarse_edge_target] = coarse_edge;
        } else {
                //we have seen this target node before and we know its postition in our
                //edge array of the graph. So we update the weight of the edge!
                EdgeWeight new_edge_weight = coarser.getEdgeWeight(edge_pos) + G.getEdgeWeight(e); 
                coarser.setEdgeWeight(edge_pos, new_edge_weight);                               
        }
}



#endif /* end of include guard: CONTRACTION_VIXZ9K0F */
