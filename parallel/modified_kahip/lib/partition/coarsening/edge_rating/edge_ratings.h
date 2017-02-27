/******************************************************************************
 * edge_ratings.h
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

#ifndef EDGE_RATING_FUNCTIONS_FUCW7H6Y
#define EDGE_RATING_FUNCTIONS_FUCW7H6Y

#include "data_structure/graph_access.h"
#include "partition_config.h"       

class edge_ratings {
public:
        edge_ratings(const PartitionConfig & partition_config);
        virtual ~edge_ratings();

        void rate(graph_access & G, unsigned level);
        void rate_expansion_star_2(graph_access & G);
        void rate_expansion_star(graph_access & G);
        void rate_expansion_star_2_algdist(graph_access & G);
        void rate_inner_outer(graph_access & G);
        void rate_pseudogeom(graph_access & G); 
        void compute_algdist(graph_access & G, std::vector<float> & dist); 
private:
        const PartitionConfig & partition_config;
};


#endif /* end of include guard: EDGE_RATING_FUNCTIONS_FUCW7H6Y */
