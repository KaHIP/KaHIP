/******************************************************************************
 * exchanger.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef EXCHANGER_YPB6QKNL
#define EXCHANGER_YPB6QKNL

#include <mpi.h>

#include "data_structure/graph_access.h"
#include "parallel_mh/population.h"
#include "partition_config.h"
#include "tools/quality_metrics.h"

class exchanger {
public:
        exchanger( MPI_Comm communicator );
        virtual ~exchanger();

        void diversify_population( PartitionConfig & config, graph_access & G, population & island, bool replace );
        void quick_start( PartitionConfig & config,  graph_access & G, population & island );
        void push_best( PartitionConfig & config,  graph_access & G, population & island );
        void recv_incoming( PartitionConfig & config,  graph_access & G, population & island );

private:
        void exchange_individum(const PartitionConfig & config, 
                                graph_access & G, 
                                int & from, 
                                int & rank, 
                                int & to, 
                                Individuum & in, Individuum & out);

        std::vector< int* >          m_partition_map_buffers;
        std::vector< MPI_Request* > m_request_pointers;
        std::vector<bool>            m_allready_send_to;

        int m_prev_best_objective;
        int m_max_num_pushes;
        int m_cur_num_pushes;

        MPI_Comm m_communicator;

        quality_metrics m_qm;
};



#endif /* end of include guard: EXCHANGER_YPB6QKNL */
