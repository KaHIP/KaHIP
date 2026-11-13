/******************************************************************************
 * kway_graph_refinement_commons.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef KWAY_GRAPH_REFINEMENT_COMMONS_PVGY97EW
#define KWAY_GRAPH_REFINEMENT_COMMONS_PVGY97EW

#include <algorithm>
#include <vector>

#include "data_structure/priority_queues/priority_queue_interface.h"
#include "definitions.h"
#include "random_functions.h"
#include "uncoarsening/refinement/refinement.h"
#include "uncoarsening/refinement/quotient_graph_refinement/2way_fm_refinement/vertex_moved_hashtable.h"

// A start node of a localized k-way search together with the gain computed
// for it just before the search and the number of random tie-breaking draws
// that computation made (the gain does not depend on the draws).
// Nodes of smaller degree are not cached: a full scan is about as cheap as a cache lookup and they
// are rarely queried more than once per round.
const EdgeWeight CACHE_MIN_DEGREE = 32;

struct kway_start_node {
        NodeID   node;
        Gain     gain;
        unsigned ties;
};

class kway_graph_refinement_commons  {
        public:

                kway_graph_refinement_commons( PartitionConfig & config );
                virtual ~kway_graph_refinement_commons();

                void init( PartitionConfig & config );

                bool incident_to_more_than_two_partitions(graph_access & G, NodeID & node);

                EdgeWeight compute_gain(graph_access & G, 
                                        NodeID & node, 
                                        PartitionID & max_gainer, 
                                        EdgeWeight & ext_degree);

                EdgeWeight compute_gain(graph_access & G, 
                                        NodeID & node, 
                                        PartitionID & max_gainer, 
                                        EdgeWeight & ext_degree,
                                        unsigned & ties);

                bool int_ext_degree( graph_access & G, 
                                     const NodeID & node,
                                     const PartitionID lhs,
                                     const PartitionID rhs,
                                     EdgeWeight & int_degree,
                                     EdgeWeight & ext_degree);

                inline unsigned getUnderlyingK();

                // Per-round connectivity cache (rule b): for each node queried during a
                // k-way refinement round, the summed edge weight towards each adjacent block
                // is kept and updated incrementally when a neighbor moves, so repeated gain
                // queries cost O(#adjacent blocks) instead of O(degree).
                inline void begin_cache_round(graph_access & G);
                inline Gain compute_gain_cached(graph_access & G,
                                                NodeID node,
                                                PartitionID & max_gainer,
                                                EdgeWeight & ext_degree);
                inline void cache_neighbor_moved(NodeID u, PartitionID from, PartitionID to, EdgeWeight w);

        private:
                struct cache_hdr {
                        unsigned stamp;
                        unsigned off;
                        unsigned size;
                        unsigned cap;
                };
                struct cache_entry {
                        PartitionID block;
                        EdgeWeight  conn;
                };
                // Shared by all instances (single-threaded): a refinement core and its commons are
                // created per localized search, so per-instance arrays of size n would cost O(n) each.
                struct cache_store {
                        std::vector<cache_hdr>   hdr;
                        std::vector<cache_entry> pool;
                        unsigned                 round = 0;
                };
                static cache_store & cache() { static cache_store s; return s; }

                //for efficient computation of internal and external degrees
                struct round_struct {
                        unsigned round;
                        EdgeWeight local_degree;
                };

                std::vector<round_struct>                    m_local_degrees;
                unsigned                                     m_round;
};

inline unsigned kway_graph_refinement_commons::getUnderlyingK() {
        return m_local_degrees.size();
}

inline void kway_graph_refinement_commons::init(PartitionConfig & config) {
        m_local_degrees.resize(config.k);
        for( PartitionID i = 0; i < config.k; i++) {
                m_local_degrees[i].round        = 0;
                m_local_degrees[i].local_degree = 0;
        }

        m_round = 0;//needed for the computation of internal and external degrees
}

inline void kway_graph_refinement_commons::begin_cache_round(graph_access & G) {
        cache_store & cs_ = cache();
        if(cs_.hdr.size() < G.number_of_nodes()) {
                cache_hdr h; h.stamp = 0; h.off = 0; h.size = 0; h.cap = 0;
                cs_.hdr.resize(G.number_of_nodes(), h);
        }
        cs_.round++;
        if(cs_.round == 0) {
                for(size_t i = 0; i < cs_.hdr.size(); i++) cs_.hdr[i].stamp = 0;
                cs_.round = 1;
        }
        cs_.pool.clear();
}

inline Gain kway_graph_refinement_commons::compute_gain_cached(graph_access & G,
                                                               NodeID node,
                                                               PartitionID & max_gainer,
                                                               EdgeWeight & ext_degree) {
        cache_store & cs_ = cache();
        if(G.getNodeDegree(node) < CACHE_MIN_DEGREE) {
                return compute_gain(G, node, max_gainer, ext_degree);
        }
        cache_hdr & h = cs_.hdr[node];
        if(h.stamp != cs_.round) {
                // build: aggregate the edge weights per adjacent block
                h.stamp = cs_.round;
                h.off   = cs_.pool.size();
                m_round++;
                forall_out_edges(G, e, node) {
                        PartitionID b = G.getPartitionIndex(G.getEdgeTarget(e));
                        EdgeWeight  w = G.getEdgeWeight(e);
                        if(m_local_degrees[b].round == m_round) {
                                cs_.pool[m_local_degrees[b].local_degree].conn += w;
                        } else {
                                m_local_degrees[b].round        = m_round;
                                m_local_degrees[b].local_degree = cs_.pool.size();
                                cache_entry ce; ce.block = b; ce.conn = w;
                                cs_.pool.push_back(ce);
                        }
                } endfor
                h.size = cs_.pool.size() - h.off;
                // the number of adjacent blocks with positive connectivity is at most min(degree, k)
                h.cap = std::max<unsigned>(h.size, std::min<size_t>(G.getNodeDegree(node), m_local_degrees.size()));
                cache_entry ce; ce.block = 0; ce.conn = 0;
                cs_.pool.resize(h.off + h.cap, ce);
        }

        PartitionID source_partition = G.getPartitionIndex(node);
        EdgeWeight max_degree = 0;
        EdgeWeight own_degree = 0;
        max_gainer = INVALID_PARTITION;
        const cache_entry * entries = &cs_.pool[0] + h.off;
        for(unsigned i = 0; i < h.size; i++) {
                const PartitionID b = entries[i].block;
                const EdgeWeight  c = entries[i].conn;
                if(b == source_partition) { own_degree = c; continue; }
                if(c > max_degree) {
                        max_degree = c;
                        max_gainer = b;
                } else if(c == max_degree) {
                        if(random_functions::nextBool()) {
                                max_gainer = b;
                        }
                }
        }
        ext_degree = max_gainer != INVALID_PARTITION ? max_degree : 0;
        return max_degree - own_degree;
}

inline void kway_graph_refinement_commons::cache_neighbor_moved(NodeID u, PartitionID from, PartitionID to, EdgeWeight w) {
        cache_store & cs_ = cache();
        cache_hdr & h = cs_.hdr[u];
        if(h.stamp != cs_.round || w == 0) return;
        cache_entry * entries = &cs_.pool[0] + h.off;
        bool found_to = false;
        for(unsigned i = 0; i < h.size; i++) {
                if(entries[i].block == from) {
                        entries[i].conn -= w;
                        if(entries[i].conn == 0) {
                                entries[i] = entries[h.size - 1];
                                h.size--;
                                i--;
                        }
                } else if(entries[i].block == to) {
                        entries[i].conn += w;
                        found_to = true;
                }
        }
        if(!found_to) {
                if(h.size == h.cap) { // only possible with zero-weight edges: relocate
                        unsigned new_off = cs_.pool.size();
                        cache_entry ce; ce.block = 0; ce.conn = 0;
                        cs_.pool.resize(new_off + 2 * h.cap + 1, ce);
                        for(unsigned i = 0; i < h.size; i++) cs_.pool[new_off + i] = cs_.pool[h.off + i];
                        h.off = new_off;
                        h.cap = 2 * h.cap + 1;
                        entries = &cs_.pool[0] + h.off;
                }
                entries[h.size].block = to;
                entries[h.size].conn  = w;
                h.size++;
        }
}

inline bool kway_graph_refinement_commons::incident_to_more_than_two_partitions(graph_access & G, NodeID & node) {
        bool ret_value = false;
        PartitionID own_partition = G.getPartitionIndex(node);
        PartitionID second_partition = INVALID_PARTITION;

        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                PartitionID target_partition = G.getPartitionIndex(target);
                if(target_partition != own_partition) {
                        if(second_partition == INVALID_PARTITION) {
                                second_partition = target_partition;
                        } else if(target_partition != second_partition) {
                                ret_value = true;
                                break;
                        }
                }

        } endfor

        return ret_value;
}

inline bool kway_graph_refinement_commons::int_ext_degree( graph_access & G, 
                                                           const NodeID & node,
                                                           const PartitionID lhs,
                                                           const PartitionID rhs,
                                                           EdgeWeight & int_degree,
                                                           EdgeWeight & ext_degree) {


        ASSERT_EQ(lhs, G.getPartitionIndex(node));

        int_degree               = 0;
        ext_degree               = 0;
        bool update_is_difficult = false;

        forall_out_edges(G, e, node) {
                NodeID target                 = G.getEdgeTarget(e);
                PartitionID targets_partition = G.getPartitionIndex(target);

                if(targets_partition == lhs) {
                        int_degree += G.getEdgeWeight(e); 
                } else if(targets_partition == rhs) {
                        ext_degree += G.getEdgeWeight(e);
                }

                if(targets_partition != lhs && targets_partition != rhs) {
                        update_is_difficult = true;
                } 
        } endfor

        return update_is_difficult;
}

inline Gain kway_graph_refinement_commons::compute_gain(graph_access & G, 
                                                        NodeID & node, 
                                                        PartitionID & max_gainer, 
                                                        EdgeWeight & ext_degree) {
        unsigned ties;
        return compute_gain(G, node, max_gainer, ext_degree, ties);
}

inline Gain kway_graph_refinement_commons::compute_gain(graph_access & G, 
                                                        NodeID & node, 
                                                        PartitionID & max_gainer, 
                                                        EdgeWeight & ext_degree,
                                                        unsigned & ties) {
        ties = 0;
        //for all incident partitions compute gain
        //return max gain and max_gainer partition
        PartitionID source_partition = G.getPartitionIndex(node);
        EdgeWeight max_degree        = 0;
        max_gainer                   = INVALID_PARTITION;

        m_round++;//can become zero again
        forall_out_edges(G, e, node) {
                NodeID target                = G.getEdgeTarget(e);
                PartitionID target_partition = G.getPartitionIndex(target);

                if(m_local_degrees[target_partition].round == m_round) {
                        m_local_degrees[target_partition].local_degree += G.getEdgeWeight(e);
                } else {
                        m_local_degrees[target_partition].local_degree = G.getEdgeWeight(e);
                        m_local_degrees[target_partition].round = m_round;
                }


                if(m_local_degrees[target_partition].local_degree >= max_degree && target_partition != source_partition) {
                        if(m_local_degrees[target_partition].local_degree > max_degree) {
                                max_degree = m_local_degrees[target_partition].local_degree;
                                max_gainer = target_partition;
                        } else {
                                //break ties randomly
                                ties++;
                                bool accept = random_functions::nextBool();
                                if(accept) {
                                        max_degree = m_local_degrees[target_partition].local_degree;
                                        max_gainer = target_partition;
                                }
                        }
                }
        } endfor

        if(max_gainer != INVALID_PARTITION) {
                ext_degree = max_degree;
        } else {
                ext_degree = 0;
        }

        if(m_local_degrees[source_partition].round != m_round) {
                m_local_degrees[source_partition].local_degree = 0;
        } 

        return max_degree-m_local_degrees[source_partition].local_degree;
}


#endif /* end of include guard: KWAY_GRAPH_REFINEMENT_COMMONS_PVGY97EW */

