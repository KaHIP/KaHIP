/******************************************************************************
 * gpa_matching.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <algorithm>
#include <deque>

#include "compare_rating.h"
#include "gpa_matching.h"
#include "macros_assertions.h"
#include "random_functions.h"

gpa_matching::gpa_matching() {

}

gpa_matching::~gpa_matching() {

}


void gpa_matching::match(const PartitionConfig & partition_config, 
                graph_access & G, 
                Matching & edge_matching, 
                CoarseMapping & coarse_mapping, 
                NodeID & no_of_coarse_vertices,
                NodePermutationMap & permutation) {
        PRINT(std::cout<< "matching using gpa" << std::endl;)
        permutation.resize(G.number_of_nodes());
        edge_matching.resize(G.number_of_nodes());
        coarse_mapping.resize(G.number_of_nodes());

        std::vector<EdgeID> edge_permutation;
        edge_permutation.reserve(G.number_of_edges());
        std::vector<NodeID> sources(G.number_of_edges());

        init(G, partition_config, permutation, edge_matching, edge_permutation, sources);

        //permutation of the edges for random tie breaking
        if(partition_config.edge_rating_tiebreaking) {
                PartitionConfig gpa_perm_config     = partition_config;
                gpa_perm_config.permutation_quality = PERMUTATION_QUALITY_GOOD;
                random_functions::permutate_entries(gpa_perm_config, edge_permutation, false);
        }

        compare_rating cmp(&G);
        std::sort(edge_permutation.begin(), edge_permutation.end(), cmp);

        path_set pathset(&G, &partition_config);

        //grow the paths
        forall_edges(G, e) {
                EdgeID curEdge = edge_permutation[e];
                NodeID source  = sources[curEdge];
                NodeID target  = G.getEdgeTarget(curEdge); 
                if(target < source) continue; // get rid of double edges

                if(G.getEdgeRating(curEdge) == 0.0) {
                        continue;
                }

                //max vertex weight constraint
                if(G.getNodeWeight(source) + G.getNodeWeight(target) > partition_config.max_vertex_weight) {
                        continue;
                }

                if( partition_config.combine ) {
                        if(G.getSecondPartitionIndex(source) != G.getSecondPartitionIndex(target) ) {
                                continue;
                        }
                }

                pathset.add_if_applicable(source, curEdge);
        } endfor 

        extract_paths_apply_matching(G, sources, edge_matching, pathset); 

        // all matched pairs are now in edge_matching 
        // now construct the coarsemapping
        no_of_coarse_vertices = 0;
        if(!partition_config.graph_allready_partitioned) {
                forall_nodes(G, n) {
                        if(partition_config.combine) {
                                if(G.getSecondPartitionIndex(n) != G.getSecondPartitionIndex(edge_matching[n])) {
                                        // v cycle... they shouldnt be contraced
                                        edge_matching[n] = n;
                                }
                        }

                        if( n < edge_matching[n]) {
                                coarse_mapping[n]                = no_of_coarse_vertices;
                                coarse_mapping[edge_matching[n]] = no_of_coarse_vertices;
                                no_of_coarse_vertices++;
                        } else if(n == edge_matching[n]) {
                                coarse_mapping[n] = no_of_coarse_vertices;
                                no_of_coarse_vertices++;
                        }

                } endfor
         } else {
                 forall_nodes(G, n) {
                        if(G.getPartitionIndex(n) != G.getPartitionIndex(edge_matching[n])) {
                                // v cycle... they shouldnt be contraced
                                edge_matching[n] = n;
                        }
                        
                        if(partition_config.combine) {
                                if(G.getSecondPartitionIndex(n) != G.getSecondPartitionIndex(edge_matching[n])) {
                                        // v cycle... they shouldnt be contraced
                                        edge_matching[n] = n;
                                }
                        }


                        if( n < edge_matching[n]) {
                                coarse_mapping[n]                = no_of_coarse_vertices;
                                coarse_mapping[edge_matching[n]] = no_of_coarse_vertices;
                                no_of_coarse_vertices++;
                        } else if(n == edge_matching[n]) {
                                coarse_mapping[n] = no_of_coarse_vertices;
                                no_of_coarse_vertices++;
                        }

                } endfor
         }
}

void gpa_matching::init(graph_access & G, 
                        const PartitionConfig & partition_config, 
                        NodePermutationMap & permutation, 
                        Matching & edge_matching, 
                        std::vector<EdgeID>  & edge_permutation, 
                        std::vector<NodeID> & sources) {

        forall_nodes(G, n) {
                permutation[n]   = n;
                edge_matching[n] = n;

                forall_out_edges(G, e, n) {
                        sources[e] = n;   
                        edge_permutation.push_back(e);

                        if(partition_config.edge_rating == WEIGHT) {
                                // in that case we need to copy it
                                G.setEdgeRating(e, G.getEdgeWeight(e));
                        }

                } endfor
        } endfor


}
void gpa_matching::extract_paths_apply_matching(graph_access & G, 
                                                std::vector<NodeID> & sources,
                                                Matching & edge_matching, 
                                                path_set & pathset) {
        // extract the paths in the path set into lists of edges.
        // then, apply the dynamic programming max weight function to them. Apply 
        // the matched edges.
        EdgeRatingType matching_rating, second_matching_rating;

        forall_nodes(G, n) {
                const path & p = pathset.get_path(n);

                if(not p.is_active()) {
                        continue;
                }
                if(p.get_tail() != n) {
                        continue;
                }
                if(p.get_length() == 0) {
                        continue;
                }

                if(p.get_head() == p.get_tail()) {
                        // ********************************
                        // handling cycles 
                        // ********************************
                        std::vector<EdgeID> a_matching, a_second_matching;
                        std::deque<EdgeID> unpacked_cycle;
                        unpack_path(p, pathset, unpacked_cycle);

                        EdgeID first = unpacked_cycle.front();
                        unpacked_cycle.pop_front();

                        maximum_weight_matching(G, 
                                        unpacked_cycle, 
                                        a_matching,
                                        matching_rating);

                        unpacked_cycle.push_front(first); 
                        EdgeID last = unpacked_cycle.back();
                        unpacked_cycle.pop_back();

                        maximum_weight_matching(G, 
                                        unpacked_cycle, 
                                        a_second_matching,
                                        second_matching_rating);

                        unpacked_cycle.push_back(last);

                        if(matching_rating > second_matching_rating) {
                                //apply first matching 
                                apply_matching(G, a_matching, sources, edge_matching); 
                        } else {
                                //apply second matching
                                apply_matching(G, a_second_matching, sources, edge_matching); 
                        }
                } else { 
                        // ********************************
                        // handling paths 
                        // ********************************
                        std::vector<EdgeID> a_matching;
                        std::vector<EdgeID> unpacked_path;

                        if(p.get_length() == 1) {
                                //match them directly
                                EdgeID e = 0;
                                if(pathset.next_vertex(p.get_tail()) == p.get_head()) {
                                        e = pathset.edge_to_next(p.get_tail());
                                } else {
                                        e = pathset.edge_to_prev(p.get_tail());
                                        ASSERT_TRUE( pathset.prev_vertex(p.get_tail()) == p.get_head() );
                                }

                                NodeID source = sources[e];
                                NodeID target = G.getEdgeTarget(e);

                                edge_matching[source] = target;
                                edge_matching[target] = source;

                                continue;                
                        }
                        unpack_path(p, pathset, unpacked_path);
                        //dump_unpacked_path(G, unpacked_path, sources); 

                        EdgeRatingType final_rating = 0;
                        maximum_weight_matching(G, unpacked_path, a_matching, final_rating);

                        //apply matched edges
                        apply_matching(G, a_matching, sources, edge_matching); 
                } 
        } endfor
}


void gpa_matching::apply_matching(graph_access & G,
                                  std::vector<EdgeID> & matched_edges, 
                                  std::vector<NodeID> & sources,
                                  Matching & edge_matching) {

        //apply matched edges 
        for( unsigned i = 0; i < matched_edges.size(); i++) {
                EdgeID e = matched_edges[i];
                NodeID source = sources[e];
                NodeID target = G.getEdgeTarget(e);

                edge_matching[source] = target;
                edge_matching[target] = source;
        }

}
template <typename VectorOrDeque> 
void gpa_matching::unpack_path(const path & p, 
                               const path_set & pathset, 
                               VectorOrDeque & unpacked_path ) {

        NodeID head = p.get_head(); 
        NodeID prev = p.get_tail();
        NodeID next;
        NodeID current = prev; 

        if(prev == head) {
                //special case: the given path is a cycle
                current = pathset.next_vertex(prev);
                unpacked_path.push_back(pathset.edge_to_next(prev));
        }

        while(current != head) {
                if(pathset.next_vertex(current) == prev) {
                        next =  pathset.prev_vertex(current);              
                        unpacked_path.push_back(pathset.edge_to_prev(current));
                } else {
                        next =  pathset.next_vertex(current);              
                        unpacked_path.push_back(pathset.edge_to_next(current));
                }
                prev = current;
                current = next;
        }
}

template <typename VectorOrDeque> 
void gpa_matching::maximum_weight_matching( graph_access & G,
                VectorOrDeque & unpacked_path, 
                std::vector<EdgeID> & matched_edges,
                EdgeRatingType & final_rating) {

        unsigned k = unpacked_path.size();
        if( k == 1 ) {
              matched_edges.push_back(unpacked_path[0]);
              return;
        }

        std::vector<EdgeRatingType> ratings(k, 0.0);
        std::vector< bool > decision(k, false);

        ratings[0] = G.getEdgeRating(unpacked_path[0]);
        ratings[1] = G.getEdgeRating(unpacked_path[1]);

        decision[0] = true;
        if(ratings[0] < ratings[1]) {
                decision[1] = true;
        }
        //build up the decision vector
        for( EdgeID i = 2; i < k; i++) {
                ASSERT_TRUE(unpacked_path[i] < G.number_of_edges());
                EdgeRatingType curRating = G.getEdgeRating(unpacked_path[i]);
                if( curRating + ratings[i-2] > ratings[i-1] ) {
                        decision[i] = true;
                        ratings[i]  = curRating + ratings[i-2];
                } else {
                        decision[i] = false;
                        ratings[i]  = ratings[i-1];
                } 
        }

        if(decision[k-1]) {
                final_rating = ratings[k-1];
        } else {
                final_rating = ratings[k-2];
        }
        //construct optimal solution 
        for(int i = k-1; i >= 0;) {
                if(decision[i]) {
                        matched_edges.push_back(unpacked_path[i]);
                        i-=2;        
                } else {
                        i-=1;
                }
        }
} 

template <typename VectorOrDeque> 
void gpa_matching::dump_unpacked_path( graph_access & G,
                VectorOrDeque & unpacked_path,
                std::vector<NodeID>& sources) {
        //dump the path
        for( unsigned i = 0; i < unpacked_path.size(); i++) {
                EdgeID e = unpacked_path[i];
                std::cout << "(" << sources[e] << " " << G.getEdgeTarget(e) << ") ";
        }
        std::cout << std::endl;


}
