/******************************************************************************
 * path_set.cpp 
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

#include "path_set.h"

path_set::path_set( graph_access * G_, const PartitionConfig * config_ ): pG(G_), config(config_),
                                         m_no_of_paths(pG->number_of_nodes()), 
                                         m_vertex_to_path(m_no_of_paths),
                                         m_paths(m_no_of_paths),
                                         m_next(m_no_of_paths), 
                                         m_prev(m_no_of_paths),
                                         m_next_edge(m_no_of_paths, UNDEFINED_EDGE),
                                         m_prev_edge(m_no_of_paths, UNDEFINED_EDGE) {
       
       graph_access & G = *pG;        
       forall_nodes(G, node) {
               m_paths[node].init(node);
               m_vertex_to_path[node] = node;
               m_next[node]           = node;
               m_prev[node]           = node;
       } endfor
       
}

path_set::~path_set() {
                
}

