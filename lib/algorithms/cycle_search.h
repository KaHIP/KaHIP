/******************************************************************************
 * cycle_search.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef CYCLE_SEARCH_IO23844C
#define CYCLE_SEARCH_IO23844C

#include "data_structure/graph_access.h"

class cycle_search {
public:
        cycle_search();
        virtual ~cycle_search();

        void find_random_cycle(graph_access & G, std::vector<NodeID> & cycle);

        //returns true if a negative cycle was found, else false
        bool find_negative_cycle(graph_access & G, NodeID & start, std::vector<NodeID> & cycle);
        
        bool find_zero_weight_cycle(graph_access & G, NodeID & start, std::vector<NodeID> & cycle); 

        bool find_shortest_path(graph_access & G, NodeID & start, NodeID & dest, std::vector<NodeID> & cycle); 

        static double total_time;
private:

        bool negative_cycle_detection(graph_access & G, 
                                      NodeID & start, 
                                      std::vector<EdgeWeight> & distance, 
                                      std::vector<NodeID> & parent, 
                                      std::vector<NodeID> & cycle);

        int bellman_ford_with_subtree_disassembly_and_updates(graph_access & G, 
                                                              NodeID & start, 
                                                              std::vector<EdgeWeight> & distance, 
                                                              std::vector<NodeID> & parent, 
                                                              std::vector<NodeID> & cycle);
};


#endif /* end of include guard: CYCLE_SEARCH_IO23844C */
