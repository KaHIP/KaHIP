/******************************************************************************
 * distributed_partitioner.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef DISTRIBUTED_PARTITIONER_ZYL2XF6R
#define DISTRIBUTED_PARTITIONER_ZYL2XF6R

#include <vector>
#include "partition_config.h"
#include "data_structure/parallel_graph_access.h"
#include "stop_rule.h"

class distributed_partitioner {
public:
        distributed_partitioner();
        virtual ~distributed_partitioner();

        void perform_partitioning( PPartitionConfig & config, parallel_graph_access & G);
        void perform_recursive_partitioning( PPartitionConfig & config, parallel_graph_access & G);

        void perform_partitioning( MPI_Comm comm, PPartitionConfig & partition_config, parallel_graph_access & G);
        void perform_recursive_partitioning( MPI_Comm comm, PPartitionConfig & partition_config, parallel_graph_access & G);

        void check( MPI_Comm comm, PPartitionConfig & config, parallel_graph_access & G);
        void check_labels( MPI_Comm comm, PPartitionConfig & config, parallel_graph_access & G);
        static void generate_random_choices( PPartitionConfig & config ) ;
private: 
        void vcycle( MPI_Comm communicator, PPartitionConfig & config, parallel_graph_access & G );

        stop_rule contraction_stop_decision;
        NodeWeight m_total_graph_weight;
        NodeID m_cur_rnd_choice;

        static std::vector< NodeID > m_cf;
        static std::vector< NodeID > m_sf;
        static std::vector< NodeID > m_lic;
        int m_level;
        int m_cycle;
        timer m_t; 
};


#endif /* end of include guard: DISTRIBUTED_PARTITIONER_ZYL2XF6R */
