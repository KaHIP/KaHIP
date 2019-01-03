/******************************************************************************
 * strongly_connected_components.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <stack>
#include <vector>

#include "strongly_connected_components.h"

strongly_connected_components::strongly_connected_components() {

}

strongly_connected_components::~strongly_connected_components() {

}

int strongly_connected_components::strong_components( graph_access & G, std::vector<int> & comp_num) {

        std::stack<NodeID> unfinished;
        std::stack<NodeID> roots;

        std::vector<int> dfsnum(G.number_of_nodes(), -1);
        m_dfscount   = 0;
        m_comp_count = 0;

        forall_nodes(G, node) {
                comp_num[node] = -1;
        } endfor

        forall_nodes(G, node) {
                if(dfsnum[node] == -1) {
                        scc_dfs(node, G, dfsnum, comp_num, unfinished, roots); 
                }
        } endfor
        return m_comp_count;
}

void strongly_connected_components::scc_dfs(NodeID node, graph_access & G, 
                                            std::vector<int> & dfsnum, 
                                            std::vector<int> & comp_num,
                                            std::stack<NodeID> & unfinished, 
                                            std::stack<NodeID> & roots){ 
        dfsnum[node] = m_dfscount++;

        //make node a tentative scc of its own
        unfinished.push(node);
        roots.push(node);

        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                //explore edge (node, target)
                if(dfsnum[target] == -1) {
                        scc_dfs(target, G, dfsnum, comp_num, unfinished, roots); 
                } else if( comp_num[target] == -1) {
                        //merge scc's
                        while( dfsnum[roots.top()] > dfsnum[target] ) roots.pop();
                }

        } endfor

        //return from call of node node
        NodeID w;
        if(node == roots.top()) {
                do {
                        w = unfinished.top(); 
                        unfinished.pop();
                        comp_num[w] = m_comp_count;
                } while( w != node );
                m_comp_count++;
                roots.pop();
        } 

}






