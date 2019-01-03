/******************************************************************************
 * gpa_matching.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef GPA_MATCHING_NXLQ0SIT
#define GPA_MATCHING_NXLQ0SIT

#include "coarsening/matching/matching.h"
#include "path.h"
#include "path_set.h"

class gpa_matching : public matching{
        public:
                gpa_matching( );
                virtual ~gpa_matching();

                void match(const PartitionConfig & config, 
                                graph_access & G, 
                                Matching & _matching, 
                                CoarseMapping & coarse_mapping, 
                                NodeID & no_of_coarse_vertices,
                                NodePermutationMap & permutation);
        private:
                void init(graph_access & G, 
                          const PartitionConfig & partition_config, 
                          NodePermutationMap & permutation, 
                          Matching & edge_matching,
                          std::vector<EdgeID>  & edge_permutation,
                          std::vector<NodeID> & sources); 

                void extract_paths_apply_matching( graph_access & G, 
                                                   std::vector<NodeID> & sources,
                                                   Matching & edge_matching, 
                                                   path_set & pathset); 

                template <typename VectorOrDeque> 
                        void unpack_path(const path & the_path, 
                                         const path_set & pathset,  
                                         VectorOrDeque & a_path);

                template <typename VectorOrDeque> 
                        void maximum_weight_matching( graph_access & G,
                                                      VectorOrDeque & unpacked_path, 
                                                      std::vector<EdgeID> & matched_edges,
                                                      EdgeRatingType & final_rating);

                void apply_matching( graph_access & G,
                                     std::vector<EdgeID> & matched_edges,
                                     std::vector<NodeID> & sources,
                                     Matching & edge_matching);


                template <typename VectorOrDeque> 
                        void dump_unpacked_path( graph_access & G,
                                                 VectorOrDeque & unpacked_path, 
                                                 std::vector<NodeID>& sources);
};


#endif /* end of include guard: GPA_MATCHING_NXLQ0SIT */
