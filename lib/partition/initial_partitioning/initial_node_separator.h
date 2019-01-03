//
// Author: Christian Schulz <christian.schulz.phone@gmail.com>
// 

#ifndef INITIAL_NODE_SEPARATOR_VR5EPEO6
#define INITIAL_NODE_SEPARATOR_VR5EPEO6

#include "data_structure/graph_access.h"
#include "partition_config.h"

class initial_node_separator {
public:
        initial_node_separator();
        virtual ~initial_node_separator();

        // method computes an initial node separator
        void compute_node_separator( const PartitionConfig & config, graph_access & G);
        NodeWeight single_run( const PartitionConfig & config, graph_access & G);
};


#endif /* end of include guard: INITIAL_NODE_SEPARATOR_VR5EPEO6 */
