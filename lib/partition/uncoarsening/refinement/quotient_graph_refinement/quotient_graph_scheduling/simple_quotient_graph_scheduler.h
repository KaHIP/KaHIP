/******************************************************************************
 * simple_quotient_graph_scheduler.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef SIMPLE_QUOTIENT_GRAPH_SCHEDULER_YG9BEBH0
#define SIMPLE_QUOTIENT_GRAPH_SCHEDULER_YG9BEBH0

#include "partition_config.h"
#include "quotient_graph_scheduling.h"

class simple_quotient_graph_scheduler : public quotient_graph_scheduling  { 
public:
        simple_quotient_graph_scheduler(PartitionConfig & config, 
                                        QuotientGraphEdges & qgraph_edges,  
                                        unsigned int account);

        virtual ~simple_quotient_graph_scheduler();

        virtual bool hasFinished();
        virtual boundary_pair & getNext();
        virtual void pushStatistics(qgraph_edge_statistics & statistic) {};

private: 
        QuotientGraphEdges m_quotient_graph_edges_pool;
};

inline bool simple_quotient_graph_scheduler::hasFinished( ) {
        return m_quotient_graph_edges_pool.empty();                
}

inline boundary_pair & simple_quotient_graph_scheduler::getNext( ) {
        boundary_pair & ret_value = m_quotient_graph_edges_pool.back();
        m_quotient_graph_edges_pool.pop_back();
        return ret_value; 
}

#endif /* end of include guard: SIMPLE_QUOTIENT_GRAPH_SCHEDULER_YG9BEBH0 */
