/******************************************************************************
 * construct_mapping.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef CONSTRUCT_MAPPING_LTW749U0
#define CONSTRUCT_MAPPING_LTW749U0

#include "data_structure/graph_access.h"
#include "data_structure/matrix/matrix.h"
#include "partition_config.h"
#include "tools/quality_metrics.h"


class construct_mapping {
        public:
                construct_mapping();
                virtual ~construct_mapping();

                void construct_initial_mapping( PartitionConfig & config, graph_access & C, matrix & D, std::vector< NodeID > & perm_rank);

                void construct_old_growing( PartitionConfig & config, graph_access & C, matrix & D, std::vector< NodeID > & perm_rank);
                void construct_old_growing_matrix( PartitionConfig & config, matrix& C, matrix & D, std::vector< NodeID > & perm_rank);
                void construct_old_growing_faster( PartitionConfig & config, graph_access & C, matrix & D, std::vector< NodeID > & perm_rank);
                void construct_identity( PartitionConfig & config, graph_access & C, matrix & D, std::vector< NodeID > & perm_rank);
                void construct_random( PartitionConfig & config, graph_access & C, matrix & D, std::vector< NodeID > & perm_rank);
                void construct_fast_hierarchy_topdown( PartitionConfig & config, graph_access & C, matrix & D, std::vector< NodeID > & perm_rank);
                void construct_fast_hierarchy_bottomup( PartitionConfig & config, graph_access & C, matrix & D, std::vector< NodeID > & perm_rank);

        private:
                int minimumNode(std::vector<int>* nodeAttribs) {
                        int minNode = -1;
                        int minValue = INT_MAX;
                        for (unsigned int i = 0; i < nodeAttribs->size(); i++) {
                                if((*nodeAttribs)[i] < minValue) {
                                        minNode = i;
                                        minValue = (*nodeAttribs)[i];
                                }
                        }
                        return minNode;
                }

                int maximumNode(std::vector<int>* nodeAttribs) {
                        int maxNode = -1;
                        int maxValue = -1;
                        for (unsigned int i = 0; i < nodeAttribs->size(); i++) {
                                if((*nodeAttribs)[i] > maxValue) {
                                        maxNode = i;
                                        maxValue = (*nodeAttribs)[i];
                                }
                        }
                        return maxNode;
                }

                quality_metrics qm;
};


#endif /* end of include guard: CONSTRUCT_MAPPING_LTW749U0 */
