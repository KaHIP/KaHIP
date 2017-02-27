/******************************************************************************
 * balance_management.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Graph Partitioning 
 ******************************************************************************
 * Copyright (C) 2017 Christian Schulz 
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
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
