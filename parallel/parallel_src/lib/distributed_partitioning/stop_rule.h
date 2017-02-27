/******************************************************************************
 * stop_rule.h
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
