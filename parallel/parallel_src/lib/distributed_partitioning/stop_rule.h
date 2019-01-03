/******************************************************************************
 * stop_rule.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef STOP_RULE_23YOZ7GX
#define STOP_RULE_23YOZ7GX

#include "data_structure/parallel_graph_access.h"
#include "partition_config.h"

class stop_rule {
public:
        stop_rule() {} ;
        virtual ~stop_rule() {};

        bool contraction_stop( PPartitionConfig & config, parallel_graph_access & finer, parallel_graph_access & coarser) {
                if( finer.number_of_global_nodes() / (double)coarser.number_of_global_nodes() < 1.1) return true;
                if( coarser.number_of_global_nodes() < config.stop_factor*config.k) return true;
                return false;       
        }
};


#endif /* end of include guard: STOP_RULE_23YOZ7GX */
