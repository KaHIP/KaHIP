/******************************************************************************
 * constraint_label_propagation.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef SIZE_CONSTRAINT_LABEL_PROPAGATION_7SVLBKKT
#define SIZE_CONSTRAINT_LABEL_PROPAGATION_7SVLBKKT

#include <unordered_map>
#include "../matching/matching.h"

struct ensemble_pair {
        PartitionID n; // number of nodes in the graph
        PartitionID lhs;
        PartitionID rhs;
};

struct compare_ensemble_pair {
        bool operator()(const ensemble_pair pair_a, const ensemble_pair pair_b) const {
                bool eq = (pair_a.lhs == pair_b.lhs && pair_a.rhs == pair_b.rhs);
                return eq;
        }
};

struct hash_ensemble_pair{
       size_t operator()(const ensemble_pair pair) const {
                return pair.lhs*pair.n + pair.rhs;
       }
};

struct data_ensemble_pair {
        NodeID mapping;

        data_ensemble_pair() {
                mapping = 0;
        }
};

typedef std::unordered_map<const ensemble_pair, 
                                data_ensemble_pair, 
                                hash_ensemble_pair, 
                                compare_ensemble_pair> hash_ensemble;


class size_constraint_label_propagation : public matching {
        public:
                size_constraint_label_propagation();
                virtual ~size_constraint_label_propagation();

                void match(const PartitionConfig & config, 
                                graph_access & G, 
                                Matching & _matching, 
                                CoarseMapping & coarse_mapping, 
                                NodeID & no_of_coarse_vertices,
                                NodePermutationMap & permutation);


                void ensemble_clusterings(const PartitionConfig & config, 
                                graph_access & G, 
                                Matching & _matching, 
                                CoarseMapping & coarse_mapping, 
                                NodeID & no_of_coarse_vertices,
                                NodePermutationMap & permutation);

                void ensemble_two_clusterings( graph_access & G,
                                std::vector<NodeID> & lhs, 
                                std::vector<NodeID> & rhs, 
                                std::vector< NodeID > & output,
                                NodeID & no_of_coarse_vertices);

                void match_internal(const PartitionConfig & config, 
                                graph_access & G, 
                                Matching & _matching, 
                                CoarseMapping & coarse_mapping, 
                                NodeID & no_of_coarse_vertices,
                                NodePermutationMap & permutation);

                void remap_cluster_ids(const PartitionConfig & partition_config, 
                                graph_access & G,
                                std::vector<NodeID> & cluster_id, 
                                NodeID & no_of_coarse_vertices,
                                bool apply_to_graph = false); 

                void create_coarsemapping(const PartitionConfig & partition_config, 
                                graph_access & G,
                                std::vector<NodeID> & cluster_id, 
                                CoarseMapping & coarse_mapping); 

                void label_propagation(const PartitionConfig & partition_config, 
                                graph_access & G,
                                const NodeWeight & block_upperbound,
                                std::vector<NodeID> & cluster_id, // output paramter
                                NodeID & number_of_blocks); // output parameter

                void label_propagation(const PartitionConfig & partition_config, 
                                graph_access & G, 
                                std::vector<NodeWeight> & cluster_id,
                                NodeID & number_of_blocks );

        private:
                // degree node order of the graph, computed once per ensemble (it depends only on G)
                std::vector<NodeID> m_degree_order;
                bool m_use_degree_order = false; 

                // G stored in m_degree_order (visiting) order, built once per ensemble: position p holds
                // node m_degree_order[p], edge targets are positions, edges in the same order per node
                void build_permuted_graph(const PartitionConfig & partition_config, graph_access & G);
                void label_propagation_permuted(const PartitionConfig & partition_config,
                                graph_access & G,
                                const NodeWeight & block_upperbound,
                                std::vector<NodeID> & cluster_id,
                                NodeID & number_of_blocks);
                std::vector<Node> m_pnodes;
                std::vector<Edge> m_pedges;
                std::vector<PartitionID> m_ppart;
                std::vector<PartitionID> m_psecond;

};


#endif /* end of include guard: SIZE_CONSTRAINT_LABEL_PROPAGATION_7SVLBKKT */
