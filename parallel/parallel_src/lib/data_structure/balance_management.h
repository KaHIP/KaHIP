/******************************************************************************
 * balance_management.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef BALANCE_MANAGEMENT_NJRUTX5K
#define BALANCE_MANAGEMENT_NJRUTX5K

#include <vector>
#include <unordered_map>
#include "definitions.h"

/* This class gives you the amount of weight that is currently in a block.
 * The information can be inaccurate but correct after a call of 
 * update_block_sizes_globally. */

class parallel_graph_access;

class balance_management {
public:
        balance_management( parallel_graph_access * G, NodeID total_num_labels);
        virtual ~balance_management();

        virtual NodeWeight getBlockSize( PartitionID block ) = 0;
        virtual void setBlockSize( PartitionID block, NodeWeight block_size ) = 0;
        virtual void update_non_contained_block_balance( PartitionID from, PartitionID to, NodeWeight node_weight) = 0;

        // init local and total block sizes
        virtual void init() = 0;
        virtual void update() = 0;

protected:
        balance_management() {};
        parallel_graph_access * m_G;
        NodeID  m_total_num_labels;

};

#endif /* end of include guard: BALANCE_MANAGEMENT_NJRUTX5K */
