/******************************************************************************
 * strongly_connected_components.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
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
