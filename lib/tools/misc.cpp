/******************************************************************************
 * misc.cpp 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
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

#include "misc.h"
#include "quality_metrics.h"

misc::misc() {

}

misc::~misc() {

}

void misc::balance_singletons(const PartitionConfig & config, graph_access & G) {
        quality_metrics qm;
        std::vector< NodeID > singletons;
        std::vector< NodeWeight > block_sizes(config.k,0);

        forall_nodes(G, node) {
                block_sizes[G.getPartitionIndex(node)] += G.getNodeWeight(node);

                if(G.getNodeDegree(node) == 0) {
                        singletons.push_back(node);
                }
        } endfor

        // use buckets?
        for( unsigned i = 0; i < singletons.size(); i++) {
                NodeWeight min = block_sizes[0];
                PartitionID p  = 0;
                for( unsigned j = 0; j < config.k; j++) {
                        if( block_sizes[j] < min ) {
                                min = block_sizes[j];
                                p   = j;
                        }
                }

                NodeID node = singletons[i];
                block_sizes[G.getPartitionIndex(node)] -= G.getNodeWeight(node);
                block_sizes[p] += G.getNodeWeight(node);
                G.setPartitionIndex(node, p);
        }
        std::cout <<  "log> balance after assigning singletons " <<  qm.balance(G)  << std::endl;
}
