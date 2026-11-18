/******************************************************************************
 * balance_management_refinement.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef BALANCE_MANAGEMENT_REFINEMENT_ZHYKQBYB
#define BALANCE_MANAGEMENT_REFINEMENT_ZHYKQBYB

#include "balance_management.h"

class parallel_graph_access;

class balance_management_refinement : public balance_management {
public:
        balance_management_refinement( parallel_graph_access * G, NodeID num_labels);
        virtual ~balance_management_refinement();

        virtual NodeWeight getBlockSize( PartitionID block );
        virtual void setBlockSize( PartitionID block, NodeWeight block_size ) ;
        virtual void update_non_contained_block_balance( PartitionID from, PartitionID to, NodeWeight node_weight) {/*noop*/};

        virtual void init();
        virtual void update();

private:
        std::vector< NodeWeight > m_total_block_weights;
        std::vector< NodeWeight > m_local_block_weights;
};


inline
void balance_management_refinement::setBlockSize( PartitionID block, NodeWeight block_size ) {
        ULONG delta = block_size - m_total_block_weights[block];
        m_local_block_weights[block] += delta;
        m_total_block_weights[block] = block_size;
}

inline
NodeWeight balance_management_refinement::getBlockSize( PartitionID block ) {
        return m_total_block_weights[block];
}

#endif /* end of include guard: BALANCE_MANAGEMENT_REFINEMENT_ZHYKQBYB */
