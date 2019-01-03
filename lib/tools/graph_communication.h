/******************************************************************************
 * graph_communication.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef GRAPH_COMMUNICATION_J5Q2P80G
#define GRAPH_COMMUNICATION_J5Q2P80G

#include "data_structure/graph_access.h"

class graph_communication {
public:
        graph_communication();
        virtual ~graph_communication();

        void broadcast_graph( graph_access & G, unsigned root);

};


#endif /* end of include guard: GRAPH_COMMUNICATION_J5Q2P80G */
