/******************************************************************************
 * size_constraint_label_propagation.h     
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/


#include <algorithm>
#include <limits>
#include <unordered_map>

#include <sstream>
#include "../edge_rating/edge_ratings.h"
#include "../matching/gpa/gpa_matching.h"
#include "data_structure/union_find.h"
#include "node_ordering.h"
#include "partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement.h"
#include "partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_commons.h"
#include "tools/quality_metrics.h"
#include "tools/random_functions.h"
#include "io/graph_io.h"

#include "size_constraint_label_propagation.h"

size_constraint_label_propagation::size_constraint_label_propagation() {
                
}

size_constraint_label_propagation::~size_constraint_label_propagation() {
                
}

void size_constraint_label_propagation::match(const PartitionConfig & partition_config, 
                                              graph_access & G, 
                                              Matching & _matching, 
                                              CoarseMapping & coarse_mapping, 
                                              NodeID & no_of_coarse_vertices,
                                              NodePermutationMap & permutation) {
        permutation.resize(G.number_of_nodes());
        coarse_mapping.resize(G.number_of_nodes());
        no_of_coarse_vertices = 0;

        if ( partition_config.ensemble_clusterings ) {
                ensemble_clusterings(partition_config, G, _matching, coarse_mapping, no_of_coarse_vertices, permutation);
        } else {
                match_internal(partition_config, G, _matching, coarse_mapping, no_of_coarse_vertices, permutation);
        }
}

void size_constraint_label_propagation::match_internal(const PartitionConfig & partition_config, 
                                              graph_access & G, 
                                              Matching & _matching, 
                                              CoarseMapping & coarse_mapping, 
                                              NodeID & no_of_coarse_vertices,
                                              NodePermutationMap & permutation) {

        std::vector<NodeWeight> cluster_id(G.number_of_nodes());
        NodeWeight block_upperbound = ceil(partition_config.upper_bound_partition/(double)partition_config.cluster_coarsening_factor);

        label_propagation( partition_config, G, block_upperbound, cluster_id, no_of_coarse_vertices);
        create_coarsemapping( partition_config, G, cluster_id, coarse_mapping);
}

void size_constraint_label_propagation::ensemble_two_clusterings( graph_access & G, 
                                                                  std::vector<NodeID> & lhs, 
                                                                  std::vector<NodeID> & rhs, 
                                                                  std::vector< NodeID > & output,
                                                                  NodeID & no_of_coarse_vertices) {

        // Each distinct pair (lhs[node], rhs[node]) gets the next id at its first occurrence in node
        // order, as with the former hash map. The nodes are bucketed by lhs (stable), so the first
        // occurrence of every pair is found with an array indexed by rhs instead of hashing.
        // output may alias lhs or rhs, so all inputs are read before output is written.
        const NodeID n = lhs.size();
        NodeID max_lhs = 0, max_rhs = 0;
        for( NodeID node = 0; node < n; node++) {
                max_lhs = std::max(max_lhs, lhs[node]);
                max_rhs = std::max(max_rhs, rhs[node]);
        }

        std::vector<NodeID> start(n == 0 ? 1 : max_lhs + 2, 0);
        for( NodeID node = 0; node < n; node++) start[lhs[node]+1]++;
        for( NodeID c = 1; c < start.size(); c++) start[c] += start[c-1];

        std::vector<NodeID> order(n);
        for( NodeID node = 0; node < n; node++) order[start[lhs[node]]++] = node;
        // start[c] is now the end of bucket c

        const NodeID unset = std::numeric_limits<NodeID>::max();
        std::vector<NodeID> first_seen(n == 0 ? 1 : max_rhs + 1, unset);
        std::vector<NodeID> first(n);
        NodeID begin = 0;
        for( NodeID c = 0; c + 1 < start.size(); c++) {
                const NodeID end = start[c];
                for( NodeID i = begin; i < end; i++) {
                        const NodeID node = order[i];
                        NodeID & f = first_seen[rhs[node]];
                        if( f == unset ) f = node;
                        first[node] = f;
                }
                for( NodeID i = begin; i < end; i++) first_seen[rhs[order[i]]] = unset;
                begin = end;
        }

        no_of_coarse_vertices = 0;
        for( NodeID node = 0; node < n; node++) {
                output[node] = first[node] == node ? no_of_coarse_vertices++ : output[first[node]];
        }
}


void size_constraint_label_propagation::ensemble_clusterings(const PartitionConfig & partition_config, 
                                                             graph_access & G, 
                                                             Matching & _matching, 
                                                             CoarseMapping & coarse_mapping, 
                                                             NodeID & no_of_coarse_vertices,
                                                             NodePermutationMap & permutation) {
        int runs = partition_config.number_of_clusterings;
        std::vector< NodeID >  cur_cluster(G.number_of_nodes(), 0);
        std::vector< NodeID >  ensemble_cluster(G.number_of_nodes(),0);

        int new_cf = partition_config.cluster_coarsening_factor;
        for( int i = 0; i < runs; i++) {
                PartitionConfig config = partition_config;
                config.cluster_coarsening_factor = new_cf;

                NodeID cur_no_blocks = 0;
                label_propagation(config, G, cur_cluster, cur_no_blocks); 

                if( i != 0 ) {
                        ensemble_two_clusterings(G, cur_cluster, ensemble_cluster, ensemble_cluster, no_of_coarse_vertices);
                } else {
                        forall_nodes(G, node) {
                                ensemble_cluster[node] = cur_cluster[node];
                        } endfor
                        
                        no_of_coarse_vertices = cur_no_blocks;
                }
                new_cf = random_functions::nextInt(10, 30);
        }

        create_coarsemapping( partition_config, G, ensemble_cluster, coarse_mapping);


}

void size_constraint_label_propagation::label_propagation(const PartitionConfig & partition_config, 
                                                         graph_access & G, 
                                                         std::vector<NodeWeight> & cluster_id, 
                                                         NodeID & no_of_blocks ) {
        NodeWeight block_upperbound = ceil(partition_config.upper_bound_partition/(double)partition_config.cluster_coarsening_factor);

        label_propagation( partition_config, G, block_upperbound, cluster_id, no_of_blocks);
}

void size_constraint_label_propagation::label_propagation(const PartitionConfig & partition_config, 
                                                         graph_access & G, 
                                                         const NodeWeight & block_upperbound,
                                                         std::vector<NodeWeight> & cluster_id,  
                                                         NodeID & no_of_blocks) {
        // in this case the _matching paramter is not used 
        // coarse_mappng stores cluster id and the mapping (it is identical)
        std::vector<PartitionID> hash_map(G.number_of_nodes(),0);
        std::vector<NodeID> permutation(G.number_of_nodes());
        std::vector<NodeWeight> cluster_sizes(G.number_of_nodes());
        cluster_id.resize(G.number_of_nodes());

        forall_nodes(G, node) {
                cluster_sizes[node] = G.getNodeWeight(node);
                cluster_id[node]    = node;
        } endfor
        
        node_ordering n_ordering;
        n_ordering.order_nodes(partition_config, G, permutation);

        // Same computation and random number sequence as before, with loop invariants hoisted and the
        // CSR arrays read through raw pointers.
        const Node * nodes = G.raw_nodes();
        const Edge * edges = G.raw_edges();
        const bool check_partition = partition_config.graph_allready_partitioned;
        const bool check_second    = partition_config.combine;
        PartitionID * hm = hash_map.data();
        NodeWeight * cs  = cluster_sizes.data();
        NodeWeight * cid = cluster_id.data();
        // The nodes are visited in random order, so the loop is bound by memory latency; prefetch the
        // node entry three nodes ahead, its edges two ahead and the cluster ids of its neighbors (at most 32) one ahead.
        const NodeID * perm = permutation.data();
        const NodeID n      = G.number_of_nodes();
        // rule (b): a node is revisited only if a neighbor changed its cluster since the node's last
        // visit (active-node label propagation); untouched neighborhoods would give the same vote.
        std::vector<char> active(n, 1);
        char * act = active.data();
        for( int j = 0; j < partition_config.label_iterations; j++) {
                unsigned int change_counter = 0;
                forall_nodes(G, i) {
                        NodeID node = perm[i];
                        if( !act[node] ) continue;
                        act[node] = 0;
                        if( i + 3 < n ) {
                                __builtin_prefetch(&nodes[perm[i+3]]);
                        }
                        if( i + 2 < n ) {
                                __builtin_prefetch(&edges[nodes[perm[i+2]].firstEdge]);
                        }
                        if( i + 1 < n ) {
                                const NodeID next = perm[i+1];
                                const EdgeID p_begin = nodes[next].firstEdge;
                                const EdgeID p_end   = std::min(nodes[next+1].firstEdge, p_begin + 32);
                                for( EdgeID e = p_begin; e < p_end; ++e) {
                                        __builtin_prefetch(&cid[edges[e].target]);
                                }
                        }
                        const EdgeID e_begin = nodes[node].firstEdge;
                        const EdgeID e_end   = nodes[node+1].firstEdge;
                        //now move the node to the cluster that is most common in the neighborhood

                        for( EdgeID e = e_begin; e < e_end; ++e) {
                                hm[cid[edges[e].target]] += edges[e].weight;
                        }

                        //second sweep for finding max and resetting array
                        PartitionID max_block = cid[node];
                        PartitionID my_block  = cid[node];
                        const NodeWeight node_weight = nodes[node].weight;

                        PartitionID max_value = 0;
                        for( EdgeID e = e_begin; e < e_end; ++e) {
                                NodeID target             = edges[e].target;
                                PartitionID cur_block     = cid[target];
                                PartitionID cur_value     = hm[cur_block];
                                if((cur_value > max_value  || (cur_value == max_value && random_functions::nextBool())) 
                                && (cs[cur_block] + node_weight < block_upperbound || cur_block == my_block) 
                                && (!check_partition || G.getPartitionIndex(node) == G.getPartitionIndex(target))
                                && (!check_second || G.getSecondPartitionIndex(node) == G.getSecondPartitionIndex(target)))
                                {
                                        max_value = cur_value;
                                        max_block = cur_block;
                                }

                                hm[cur_block] = 0;
                        }

                        cs[cid[node]]  -= node_weight;
                        cs[max_block]  += node_weight;
                        if( cid[node] != max_block ) {
                                change_counter++;
                                for( EdgeID e = e_begin; e < e_end; ++e) act[edges[e].target] = 1;
                        }
                        cid[node]       = max_block;
                } endfor
                if( change_counter == 0 ) break;
        }

        remap_cluster_ids( partition_config, G, cluster_id, no_of_blocks);
}



void size_constraint_label_propagation::create_coarsemapping(const PartitionConfig & partition_config, 
                                                             graph_access & G,
                                                             std::vector<NodeWeight> & cluster_id,
                                                             CoarseMapping & coarse_mapping) {
        forall_nodes(G, node) {
                coarse_mapping[node] = cluster_id[node];
        } endfor
}

void size_constraint_label_propagation::remap_cluster_ids(const PartitionConfig & partition_config, 
                                                          graph_access & G,
                                                          std::vector<NodeWeight> & cluster_id,
                                                          NodeID & no_of_coarse_vertices, bool apply_to_graph) {

        PartitionID cur_no_clusters = 0;
        std::unordered_map<PartitionID, PartitionID> remap;
        forall_nodes(G, node) {
                PartitionID cur_cluster = cluster_id[node];
                //check wether we already had that
                if( remap.find( cur_cluster ) == remap.end() ) {
                        remap[cur_cluster] = cur_no_clusters++;
                }

                cluster_id[node] = remap[cur_cluster];
        } endfor

        if( apply_to_graph ) {
                forall_nodes(G, node) {
                        G.setPartitionIndex(node, cluster_id[node]);
                } endfor
                G.set_partition_count(cur_no_clusters);
        }

        no_of_coarse_vertices = cur_no_clusters;
}

