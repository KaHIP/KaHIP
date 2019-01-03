/******************************************************************************
 * graph_hierarchy.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef GRAPH_HIERACHY_UMHG74CO
#define GRAPH_HIERACHY_UMHG74CO

#include <stack>

#include "graph_access.h"
#include "uncoarsening/refinement/quotient_graph_refinement/partial_boundary.h"

class graph_hierarchy {
public:
        graph_hierarchy( );
        virtual ~graph_hierarchy();

        void push_back(graph_access * G, CoarseMapping * coarse_mapping);
        
        graph_access  * pop_finer_and_project();
        graph_access  * pop_finer_and_project_ns( PartialBoundary & separator );
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
