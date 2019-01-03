/******************************************************************************
 * tabu_search.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef TABU_SEARCH_RC6W8GX
#define TABU_SEARCH_RC6W8GX

#include "data_structure/matrix/matrix.h"
#include "data_structure/matrix/normal_matrix.h"
#include "definitions.h"
#include "uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_commons.h"
#include "uncoarsening/refinement/refinement.h"

class tabu_search : public refinement {
        public:
                tabu_search();
                virtual ~tabu_search();

                virtual EdgeWeight perform_refinement(PartitionConfig & config, 
                                                      graph_access & G, 
                                                      complete_boundary & boundary); 

	private:
		unsigned compute_tenure(unsigned iteration, unsigned max_iteration) {
			 std::vector< double > b(15,0);
			 b[0]  = 1/8.0;
			 b[1]  = 2/8.0;
			 b[2]  = 1/8.0;
			 b[3]  = 4/8.0;
			 b[4]  = 1/8.0;
			 b[5]  = 2/8.0;
			 b[6]  = 1/8.0;
			 b[7]  = 8/8.0;
			 b[8]  = 1/8.0;
			 b[9]  = 2/8.0;
			 b[10] = 1/8.0;
			 b[11] = 4/8.0;
			 b[12] = 1/8.0;
			 b[13] = 2/8.0;
			 b[14] = 1/8.0;

			 //compute i
			 unsigned i = 1;
			 unsigned x = 4*max_iteration*b[0];
			 while( true ) {
				 if( iteration >= x ) {
					x = x + 4*max_iteration*b[i%15];
					i++;
				 } else {
				 	i--;
					break;
				 }
			 }
			 return max_iteration*b[i%15];
		
		}

                kway_graph_refinement_commons* commons;
                matrix* m;
};


#endif /* end of include guard: TABU_SEARCH_RC6W8GGX */
