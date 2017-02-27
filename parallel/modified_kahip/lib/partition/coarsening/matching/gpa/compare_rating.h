/******************************************************************************
 * compare_rating.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Graph Partitioning 
 ******************************************************************************
 * Copyright (C) 2017 Christian Schulz 
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
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
