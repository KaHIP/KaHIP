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

        m_active_nodes = true;
        if ( partition_config.ensemble_clusterings ) {
                ensemble_clusterings(partition_config, G, _matching, coarse_mapping, no_of_coarse_vertices, permutation);
        } else {
                match_internal(partition_config, G, _matching, coarse_mapping, no_of_coarse_vertices, permutation);
        }
        m_active_nodes = false;
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

        // The degree order is deterministic and the graph does not change between the runs, so it is
        // sorted once (the random order draws randoms and is still computed per run).
        if( partition_config.node_ordering == DEGREE_NODEORDERING ) {
                m_degree_order.resize(G.number_of_nodes());
                node_ordering n_ordering;
                n_ordering.order_nodes(partition_config, G, m_degree_order);
                m_use_degree_order = true;
                build_permuted_graph(partition_config, G);
        }

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

        m_use_degree_order = false;
        std::vector<Node>().swap(m_pnodes);
        std::vector<Edge>().swap(m_pedges);
        std::vector<PartitionID>().swap(m_ppart);
        std::vector<PartitionID>().swap(m_psecond);
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
        if( m_use_degree_order ) { // ensemble: order and permuted copy of G built once for all clusterings
                label_propagation_permuted( partition_config, G, block_upperbound, cluster_id, no_of_blocks);
                return;
        }

        // a single clustering also runs on a copy of G in its visiting order
        reuse_store & r = reuse();
        if( r.graph == &G && partition_config.node_ordering == DEGREE_NODEORDERING
            && !partition_config.graph_allready_partitioned && !partition_config.combine ) {
                if( !r.valid ) {
                        std::vector<NodeID> permutation(G.number_of_nodes());
                        node_ordering n_ordering;
                        n_ordering.order_nodes(partition_config, G, permutation);
                        r.order.swap(permutation);
                        m_degree_order.swap(r.order);
                        build_permuted_graph(partition_config, G);
                        m_degree_order.swap(r.order);
                        r.pnodes.swap(m_pnodes);
                        r.pedges.swap(m_pedges);
                        r.valid = true;
                }
                m_degree_order.swap(r.order);
                m_pnodes.swap(r.pnodes);
                m_pedges.swap(r.pedges);
                label_propagation_permuted( partition_config, G, block_upperbound, cluster_id, no_of_blocks);
                m_degree_order.swap(r.order);
                m_pnodes.swap(r.pnodes);
                m_pedges.swap(r.pedges);
                return;
        }
        std::vector<NodeID> permutation(G.number_of_nodes());
        node_ordering n_ordering;
        n_ordering.order_nodes(partition_config, G, permutation);
        m_degree_order.swap(permutation);
        build_permuted_graph(partition_config, G);
        label_propagation_permuted( partition_config, G, block_upperbound, cluster_id, no_of_blocks);
        m_degree_order.swap(permutation);
}



void size_constraint_label_propagation::set_reuse_graph(const graph_access * G) {
        reuse_store & r = reuse();
        r.graph = G;
        r.valid = false;
        std::vector<NodeID>().swap(r.order);
        std::vector<Node>().swap(r.pnodes);
        std::vector<Edge>().swap(r.pedges);
}

void size_constraint_label_propagation::build_permuted_graph(const PartitionConfig & partition_config, graph_access & G) {
        const NodeID n = G.number_of_nodes();
        const Node * nodes = G.raw_nodes();
        const Edge * edges = G.raw_edges();
        std::vector<NodeID> pos(n);
        for( NodeID p = 0; p < n; p++) pos[m_degree_order[p]] = p;
        m_pnodes.resize(n+1);
        m_pedges.resize(G.number_of_edges());
        EdgeID cur = 0;
        for( NodeID p = 0; p < n; p++) {
                const NodeID node = m_degree_order[p];
                m_pnodes[p].firstEdge = cur;
                m_pnodes[p].weight    = nodes[node].weight;
                for( EdgeID e = nodes[node].firstEdge; e < nodes[node+1].firstEdge; ++e, ++cur) {
                        m_pedges[cur].target = pos[edges[e].target];
                        m_pedges[cur].weight = edges[e].weight;
                }
        }
        m_pnodes[n].firstEdge = cur;
        m_pnodes[n].weight    = 0;
        if( partition_config.graph_allready_partitioned ) {
                m_ppart.resize(n);
                for( NodeID p = 0; p < n; p++) m_ppart[p] = G.getPartitionIndex(m_degree_order[p]);
        }
        if( partition_config.combine ) {
                m_psecond.resize(n);
                for( NodeID p = 0; p < n; p++) m_psecond[p] = G.getSecondPartitionIndex(m_degree_order[p]);
        }
}

// The label propagation above on the permuted copy of G: position p (node m_degree_order[p]) is
// visited p-th, so the node and edge arrays are read sequentially. The clusters are labeled by
// positions instead of node ids; the computation only compares labels for equality, so it makes the
// same choices with the same random numbers, and remap_cluster_ids (first occurrence in node order)
// gives the same ids.
void size_constraint_label_propagation::label_propagation_permuted(const PartitionConfig & partition_config, 
                                                         graph_access & G, 
                                                         const NodeWeight & block_upperbound,
                                                         std::vector<NodeWeight> & cluster_id,  
                                                         NodeID & no_of_blocks) {
        struct cluster_slot {
                PartitionID votes;
                NodeWeight  size;
        };
        const NodeID n = G.number_of_nodes();
        const Node * nodes = m_pnodes.data();
        const Edge * edges = m_pedges.data();
        std::vector<cluster_slot> slots(n);
        std::vector<NodeID> pcid(n);
        for( NodeID p = 0; p < n; p++) {
                slots[p].votes = 0;
                slots[p].size  = nodes[p].weight;
                pcid[p]        = p;
        }

        const bool check_partition = partition_config.graph_allready_partitioned;
        const bool check_second    = partition_config.combine;
        const PartitionID * ppart   = m_ppart.data();
        const PartitionID * psecond = m_psecond.data();
        cluster_slot * cl = slots.data();
        NodeID * cid = pcid.data();
        const bool active_nodes = m_active_nodes;
        std::vector<char> active(n, 1);
        char * act = active.data();
        for( int j = 0; j < partition_config.label_iterations; j++) {
                unsigned int change_counter = 0;
                for( NodeID node = 0; node < n; node++) {
                        if( active_nodes ) {
                                if( !act[node] ) continue;
                                act[node] = 0;
                        }
                        if( node + 1 < n ) {
                                const EdgeID p_begin = nodes[node+1].firstEdge;
                                const EdgeID p_end   = std::min(nodes[node+2].firstEdge, p_begin + 32);
                                for( EdgeID e = p_begin; e < p_end; ++e) {
                                        __builtin_prefetch(&cid[edges[e].target]);
                                }
                        }
                        const EdgeID e_begin = nodes[node].firstEdge;
                        const EdgeID e_end   = nodes[node+1].firstEdge;

                        for( EdgeID e = e_begin; e < e_end; ++e) {
                                cl[cid[edges[e].target]].votes += edges[e].weight;
                        }

                        PartitionID max_block = cid[node];
                        PartitionID my_block  = cid[node];
                        const NodeWeight node_weight = nodes[node].weight;

                        PartitionID max_value = 0;
                        for( EdgeID e = e_begin; e < e_end; ++e) {
                                NodeID target             = edges[e].target;
                                PartitionID cur_block     = cid[target];
                                PartitionID cur_value     = cl[cur_block].votes;
                                if((cur_value > max_value  || (cur_value == max_value && random_functions::nextBool())) 
                                && (cl[cur_block].size + node_weight < block_upperbound || cur_block == my_block) 
                                && (!check_partition || ppart[node] == ppart[target])
                                && (!check_second || psecond[node] == psecond[target]))
                                {
                                        max_value = cur_value;
                                        max_block = cur_block;
                                }

                                cl[cur_block].votes = 0;
                        }

                        cl[cid[node]].size  -= node_weight;
                        cl[max_block].size  += node_weight;
                        if( cid[node] != max_block ) {
                                change_counter++;
                                if( active_nodes ) {
                                        for( EdgeID e = e_begin; e < e_end; ++e) act[edges[e].target] = 1;
                                }
                        }
                        cid[node]       = max_block;
                }
                if( active_nodes && change_counter == 0 ) break;
        }

        cluster_id.resize(n);
        for( NodeID p = 0; p < n; p++) cluster_id[m_degree_order[p]] = pcid[p];
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

        // ids in order of first occurrence, as with the former hash map; an array indexed by the old id
        PartitionID cur_no_clusters = 0;
        PartitionID max_id = 0;
        forall_nodes(G, node) {
                max_id = std::max(max_id, (PartitionID)cluster_id[node]);
        } endfor
        const PartitionID unset = std::numeric_limits<PartitionID>::max();
        std::vector<PartitionID> remap(G.number_of_nodes() == 0 ? 1 : max_id + 1, unset);
        forall_nodes(G, node) {
                PartitionID & r = remap[cluster_id[node]];
                if( r == unset ) r = cur_no_clusters++;
                cluster_id[node] = r;
        } endfor

        if( apply_to_graph ) {
                forall_nodes(G, node) {
                        G.setPartitionIndex(node, cluster_id[node]);
                } endfor
                G.set_partition_count(cur_no_clusters);
        }

        no_of_coarse_vertices = cur_no_clusters;
}

