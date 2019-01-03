/******************************************************************************
 * path_set.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
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

