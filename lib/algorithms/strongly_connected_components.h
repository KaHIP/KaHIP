/******************************************************************************
 * strongly_connected_components.h 
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

#ifndef STRONGLY_CONNECTED_COMPONENTS_7ZJ8233R
#define STRONGLY_CONNECTED_COMPONENTS_7ZJ8233R

#include <stack>
#include <vector>

#include "data_structure/graph_access.h"
#include "definitions.h"

class strongly_connected_components {
public:
        strongly_connected_components();
        virtual ~strongly_connected_components();

        int strong_components( graph_access & G, std::vector<int> & comp_num);       
        void explicit_scc_dfs(NodeID node, graph_access & G); 

private:
        int m_dfscount; 
        int m_comp_count;

        std::vector<int> m_dfsnum;
        std::vector<int> m_comp_num;
        std::stack<NodeID> m_unfinished;
        std::stack<NodeID> m_roots;
        std::stack< std::pair<NodeID,EdgeID> > iteration_stack;
};


#endif /* end of include guard: STRONGLY_CONNECTED_COMPONENTS_7ZJ8233R */
