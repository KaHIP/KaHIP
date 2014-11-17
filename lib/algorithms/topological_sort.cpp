/******************************************************************************
 * topological_sort.cpp 
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

#include <algorithm>

#include "random_functions.h"
#include "topological_sort.h"

topological_sort::topological_sort() {
                
}

topological_sort::~topological_sort() {
                
}

void topological_sort::sort( graph_access & G, std::vector<NodeID> & sorted_sequence) {
        std::vector<int> dfsnum(G.number_of_nodes(), -1);
        int dfscount = 0;

        std::vector<NodeID> nodes(G.number_of_nodes());
        random_functions::permutate_vector_good(nodes, true);

        forall_nodes(G, node) {
                NodeID curNode = nodes[node];
                if(dfsnum[curNode] == -1) {
                        sort_dfs(curNode, G, dfsnum, dfscount, sorted_sequence); 
                }
        } endfor

        std::reverse(sorted_sequence.begin(), sorted_sequence.end());
}


void topological_sort::sort_dfs(NodeID node, graph_access & G, 
                                std::vector<int> & dfsnum, 
                                int & dfscount,
                                std::vector<NodeID> & sorted_sequence){ 

        dfsnum[node] = dfscount++;
        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                //explore edge (node, target)
                if(dfsnum[target] == -1) {
                        sort_dfs(target, G, dfsnum, dfscount, sorted_sequence); 
                } 
        } endfor
       
        //return from call of node node
        sorted_sequence.push_back(node);
}
