/******************************************************************************
 * parallel_block_down_propagation.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef PARALLEL_BLOCK_DOWN_PROPAGATION_SRTCMH8F
#define PARALLEL_BLOCK_DOWN_PROPAGATION_SRTCMH8F

#include "data_structure/parallel_graph_access.h"
#include "partition_config.h"

class parallel_block_down_propagation {
public:
        parallel_block_down_propagation();
        virtual ~parallel_block_down_propagation();

        void propagate_block_down( MPI_Comm communicator, PPartitionConfig & config, 
                                   parallel_graph_access & G, 
                                   parallel_graph_access & Q);

private:

        void update_ghost_nodes_blocks( MPI_Comm communicator, parallel_graph_access & G ); 

        std::vector< std::vector< NodeID > > m_messages;
        std::vector< std::vector< NodeID > > m_send_buffers; // buffers to send messages
};


#endif /* end of include guard: PARALLEL_BLOCK_DOWN_PROPAGATION_SRTCMH8F */
