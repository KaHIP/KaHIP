/******************************************************************************
 * edge_ratings.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef EDGE_RATING_FUNCTIONS_FUCW7H6Y
#define EDGE_RATING_FUNCTIONS_FUCW7H6Y

#include "data_structure/graph_access.h"
#include "partition_config.h"       

class edge_ratings {
public:
        edge_ratings(const PartitionConfig & partition_config);
        virtual ~edge_ratings();

        void rate(graph_access & G, unsigned level);
        void rate_expansion_star_2(graph_access & G);
        void rate_expansion_star(graph_access & G);
        void rate_expansion_star_2_algdist(graph_access & G);
        void rate_inner_outer(graph_access & G);
        void rate_pseudogeom(graph_access & G); 
        void compute_algdist(graph_access & G, std::vector<float> & dist); 
        void rate_separator_addx(graph_access & G);
        void rate_separator_multx(graph_access & G);
        void rate_separator_max(graph_access & G);
        void rate_separator_log(graph_access & G);
        void rate_separator_r1(graph_access & G);
        void rate_separator_r2(graph_access & G);
        void rate_separator_r3(graph_access & G);
        void rate_separator_r4(graph_access & G);
        void rate_separator_r5(graph_access & G);
        void rate_separator_r6(graph_access & G);
        void rate_separator_r7(graph_access & G);
        void rate_separator_r8(graph_access & G);
        void rate_realweight(graph_access & G);

private:
        const PartitionConfig & partition_config;
};


#endif /* end of include guard: EDGE_RATING_FUNCTIONS_FUCW7H6Y */
