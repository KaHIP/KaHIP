/******************************************************************************
 * topological_sort.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
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
