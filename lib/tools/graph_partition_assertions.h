/******************************************************************************
 * graph_partition_assertions.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef GRAPH_PARTITION_ASSERTIONS_609QZZDM
#define GRAPH_PARTITION_ASSERTIONS_609QZZDM

#include "data_structure/graph_access.h"
#include "partition_config.h"

class graph_partition_assertions {
        public:
                graph_partition_assertions( ) {};
                virtual ~graph_partition_assertions() {};

                static bool assert_graph_has_kway_partition(const PartitionConfig & config, graph_access & G) {
                        bool* allpartsthere = new bool[config.k];
                        for(unsigned int i = 0; i < config.k; i++) {
                                allpartsthere[i] = false;
                        }

                        forall_nodes(G, n) {
                                allpartsthere[G.getPartitionIndex(n)] = true; 
                        } endfor

                        for(unsigned int i = 0; i < config.k; i++) {
                                ASSERT_TRUE(allpartsthere[i]);
                        }

                        delete[] allpartsthere;
                        return true;
                };

};


#endif /* end of include guard: GRAPH_PARTITION_ASSERTIONS_609QZZDM */
