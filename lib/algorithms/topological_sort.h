/******************************************************************************
 * topological_sort.h 
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

#ifndef TOPOLOGICAL_SORT_GB9FC2CZ
#define TOPOLOGICAL_SORT_GB9FC2CZ

#include <stack>
#include <vector>

#include "data_structure/graph_access.h"
#include "definitions.h"

class topological_sort {
public:
        topological_sort();
        virtual ~topological_sort();

        void sort( graph_access & SG, std::vector<NodeID> & sorted_sequence);

        void sort_dfs(NodeID node, graph_access & G, 
                      std::vector<int>    & dfsnum, 
                      int                 & dfscount,
                      std::vector<NodeID> & sorted_sequence);
};


#endif /* end of include guard: TOPOLOGICAL_SORT_GB9FC2CZ */
