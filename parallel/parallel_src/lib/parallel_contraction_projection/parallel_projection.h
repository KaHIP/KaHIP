/******************************************************************************
 * parallel_projection.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef PARALLEL_PROJECTION_HBRCPQ0P
#define PARALLEL_PROJECTION_HBRCPQ0P

#include "data_structure/parallel_graph_access.h"

class parallel_projection {
public:
        parallel_projection();
        virtual ~parallel_projection();

        void parallel_project( MPI_Comm communicator, parallel_graph_access & finer, parallel_graph_access & coarser );

        //initial assignment after initial partitioning
        void initial_assignment( parallel_graph_access & G, complete_graph_access & Q);
private:
        std::vector< std::vector< NodeID > > m_messages;
};




#endif /* end of include guard: PARALLEL_PROJECTION_HBRCPQ0P */
