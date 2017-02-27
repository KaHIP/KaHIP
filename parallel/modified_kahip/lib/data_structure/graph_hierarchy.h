/******************************************************************************
 * graph_hierarchy.h
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

#ifndef GRAPH_HIERACHY_UMHG74CO
#define GRAPH_HIERACHY_UMHG74CO

#include <stack>

#include "graph_access.h"

class graph_hierarchy {
public:
        graph_hierarchy( );
        virtual ~graph_hierarchy();

        void push_back(graph_access * G, CoarseMapping * coarse_mapping);
        
        graph_access  * pop_finer_and_project();
        graph_access  * get_coarsest();
        CoarseMapping * get_mapping_of_current_finer();
               
        bool isEmpty();
        unsigned int size();
private:
        //private functions
        graph_access * pop_coarsest();

        std::stack<graph_access*>   m_the_graph_hierarchy;
        std::stack<CoarseMapping*>  m_the_mappings;
        std::vector<CoarseMapping*> m_to_delete_mappings;
        std::vector<graph_access*>  m_to_delete_hierachies;
        graph_access  * m_current_coarser_graph;
        graph_access  * m_coarsest_graph;
        CoarseMapping * m_current_coarse_mapping;
};


#endif /* end of include guard: GRAPH_HIERACHY_UMHG74CO */
