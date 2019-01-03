/******************************************************************************
 * balance_management_coarsening.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef BALANCE_MANAGEMENT_COARSENING_TS6EZN5A
#define BALANCE_MANAGEMENT_COARSENING_TS6EZN5A

#include "balance_management.h"

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

private:
        std::unordered_map< PartitionID, long > m_fuzzy_block_weights;
};


inline 
NodeWeight balance_management_coarsening::getBlockSize( PartitionID block ) {
        return m_fuzzy_block_weights[block];
}

// this function is only called for local nodes
inline 
void balance_management_coarsening::setBlockSize( PartitionID block, NodeWeight block_size ) {
        if( block_size == 0 ) {
                m_fuzzy_block_weights.erase(block);
        } else {
                m_fuzzy_block_weights[block] = block_size;
        }
}

inline
void balance_management_coarsening::update_non_contained_block_balance( PartitionID from, PartitionID to, NodeWeight node_weight) {
        if( m_fuzzy_block_weights[from] == (long)node_weight) {
                m_fuzzy_block_weights.erase(from);
        } else {
                m_fuzzy_block_weights[from] -= node_weight;
        }
        if( m_fuzzy_block_weights.find( to ) == m_fuzzy_block_weights.end() ) {
                m_fuzzy_block_weights[to] = node_weight;
        } else {
                m_fuzzy_block_weights[to] += node_weight;
        }
}

#endif /* end of include guard: BALANCE_MANAGEMENT_COARSENING_TS6EZN5A */
