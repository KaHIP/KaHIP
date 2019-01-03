/******************************************************************************
 * compare_rating.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef COMPARE_RATING_750FUZ7Z
#define COMPARE_RATING_750FUZ7Z

#include "data_structure/graph_access.h"
#include "definitions.h"

class compare_rating : public std::binary_function<EdgeRatingType, EdgeRatingType, bool> {
        public:
                compare_rating(graph_access * pG) : G(pG) {};
                virtual ~compare_rating() {};

                bool operator() (const EdgeRatingType left, const EdgeRatingType right ) {
                        return G->getEdgeRating(left) > G->getEdgeRating(right);
                }

        private:
                graph_access * G;
};


#endif /* end of include guard: COMPARE_RATING_750FUZ7Z */
