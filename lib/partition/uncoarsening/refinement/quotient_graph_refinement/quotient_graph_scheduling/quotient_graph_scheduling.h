/******************************************************************************
 * quotient_graph_scheduling.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef QUOTIENT_GRAPH_SCHEDULING_NEFT9H3J
#define QUOTIENT_GRAPH_SCHEDULING_NEFT9H3J

#include "definitions.h"
#include "uncoarsening/refinement/quotient_graph_refinement/complete_boundary.h"

struct qgraph_edge_statistics {

        EdgeWeight improvement;
        bool something_changed;
        boundary_pair* pair;

        qgraph_edge_statistics(EdgeWeight _improvement, 
                               boundary_pair* bp, 
                               bool change) : improvement(_improvement), something_changed(change), pair(bp){
        }
};

class quotient_graph_scheduling {
        public:
                quotient_graph_scheduling();
                virtual ~quotient_graph_scheduling();

                virtual bool hasFinished() = 0;
                virtual boundary_pair & getNext() = 0;
                virtual void pushStatistics(qgraph_edge_statistics & statistic) = 0;

};


#endif /* end of include guard: QUOTIENT_GRAPH_SCHEDULING_NEFT9H3J */
