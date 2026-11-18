/******************************************************************************
 * balance_management_coarsening.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef BALANCE_MANAGEMENT_COARSENING_TS6EZN5A
#define BALANCE_MANAGEMENT_COARSENING_TS6EZN5A

#include <vector>

#include "balance_management.h"
#include "block_weight_map.h"

class parallel_graph_access;

class balance_management_coarsening : public balance_management {
public:
        balance_management_coarsening( parallel_graph_access * G, PartitionID num_labels );
        virtual ~balance_management_coarsening();

        virtual NodeWeight getBlockSize( PartitionID block );
        virtual void setBlockSize( PartitionID block, NodeWeight block_size );
        virtual void update_non_contained_block_balance( PartitionID from, PartitionID to, NodeWeight node_weight);

        virtual void init();
        virtual void update();

        // the weights themselves, so that parallel_graph_access can read and update them without a virtual call
        block_weight_map * weights() { return &m_fuzzy_block_weights; }

private:
        block_weight_map m_fuzzy_block_weights;
};


inline 
NodeWeight balance_management_coarsening::getBlockSize( PartitionID block ) {
        return m_fuzzy_block_weights.get(block);
}

// this function is only called for local nodes
inline 
void balance_management_coarsening::setBlockSize( PartitionID block, NodeWeight block_size ) {
        m_fuzzy_block_weights.ref(block) = block_size;
}

inline
void balance_management_coarsening::update_non_contained_block_balance( PartitionID from, PartitionID to, NodeWeight node_weight) {
        if( from == to ) return; // subtracting and adding the same weight changes nothing
        m_fuzzy_block_weights.ref(from) -= node_weight;
        m_fuzzy_block_weights.ref(to)   += node_weight;
}

#endif /* end of include guard: BALANCE_MANAGEMENT_COARSENING_TS6EZN5A */
