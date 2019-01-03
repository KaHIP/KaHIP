/******************************************************************************
 * kway_graph_refinement_commons.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef KWAY_GRAPH_REFINEMENT_COMMONS_PVGY97EW
#define KWAY_GRAPH_REFINEMENT_COMMONS_PVGY97EW

#include <vector>

#include "data_structure/priority_queues/priority_queue_interface.h"
#include "definitions.h"
#include "random_functions.h"
#include "uncoarsening/refinement/refinement.h"
#include "uncoarsening/refinement/quotient_graph_refinement/2way_fm_refinement/vertex_moved_hashtable.h"

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

                bool int_ext_degree( graph_access & G, 
                                     const NodeID & node,
                                     const PartitionID lhs,
                                     const PartitionID rhs,
                                     EdgeWeight & int_degree,
                                     EdgeWeight & ext_degree);

                inline unsigned getUnderlyingK();

        private:

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

