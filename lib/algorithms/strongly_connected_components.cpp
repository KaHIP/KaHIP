/******************************************************************************
 * strongly_connected_components.cpp 
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

#include <stack>
#include <vector>

#include "strongly_connected_components.h"

strongly_connected_components::strongly_connected_components() {

}

strongly_connected_components::~strongly_connected_components() {

}

int strongly_connected_components::strong_components( graph_access & G, std::vector<int> & comp_num) {

        m_dfsnum.resize(G.number_of_nodes());
        m_comp_num.resize(G.number_of_nodes());
        m_dfscount   = 0;
        m_comp_count = 0;

        forall_nodes(G, node) {
                //comp_num[node] = -1;
                m_comp_num[node] = -1;
                m_dfsnum[node]   = -1;
        } endfor

        forall_nodes(G, node) {
                if(m_dfsnum[node] == -1) {
                        explicit_scc_dfs(node, G); 
                }
        } endfor

        forall_nodes(G, node) {
                comp_num[node] = m_comp_num[node];
        } endfor
        
        return m_comp_count;
}

void strongly_connected_components::explicit_scc_dfs(NodeID node, graph_access & G){ 

        iteration_stack.push( std::pair<NodeID,EdgeID>( node, G.get_first_edge(node) ) );

        //make node a tentative scc of its own
        m_dfsnum[node] = m_dfscount++;
        m_unfinished.push(node);
        m_roots.push(node);

        while( !iteration_stack.empty() ) {
                NodeID current_node = iteration_stack.top().first;
                EdgeID current_edge = iteration_stack.top().second;
                iteration_stack.pop();

                forall_out_edges_starting_at(G, e, current_node, current_edge) {
                        NodeID target = G.getEdgeTarget(e);
                        //explore edge (node, target)
                        if(m_dfsnum[target] == -1) {

                                iteration_stack.push( std::pair<NodeID,EdgeID>( current_node, e ) );
                                iteration_stack.push( std::pair<NodeID,EdgeID>( target, G.get_first_edge(target) ) );

                                m_dfsnum[target] = m_dfscount++;
                                m_unfinished.push(target);
                                m_roots.push(target);
                                break;
                        } else if( m_comp_num[target] == -1) {
                                //merge scc's
                                while( m_dfsnum[m_roots.top()] > m_dfsnum[target] ) m_roots.pop();
                        }

                } endfor

                //return from call of node node
                if(current_node == m_roots.top()) {
                        NodeID w = 0;
                        do {
                                w = m_unfinished.top(); 
                                m_unfinished.pop();
                                m_comp_num[w] = m_comp_count;
                        } while( w != current_node );
                        m_comp_count++;
                        m_roots.pop();
                } 
        }

}


