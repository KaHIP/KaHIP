//
// Author: Christian Schulz <christian.schulz@kit.edu>
// 

#ifndef GREEDY_NS_LOCAL_SEARCH_P9KLE4NH
#define GREEDY_NS_LOCAL_SEARCH_P9KLE4NH

#include "definitions.h"
#include "partition_config.h"
#include "data_structure/graph_access.h"

class greedy_ns_local_search {
public:
        greedy_ns_local_search();
        virtual ~greedy_ns_local_search();

        EdgeWeight perform_refinement(const PartitionConfig & config, graph_access & G);

};


#endif /* end of include guard: GREEDY_NS_LOCAL_SEARCH_P9KLE4NH */
