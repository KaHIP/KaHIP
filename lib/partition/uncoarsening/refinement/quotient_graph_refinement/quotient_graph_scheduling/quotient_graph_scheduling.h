/******************************************************************************
 * quotient_graph_scheduling.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
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

#ifndef QUOTIENT_GRAPH_SCHEDULING_NEFT9H3J
#define QUOTIENT_GRAPH_SCHEDULING_NEFT9H3J

#include "definitions.h"
#include "uncoarsening/refinement/quotient_graph_refinement/complete_boundary.h"

struct qgraph_edge_statistics {

        EdgeWeight improvement;
        bool something_changed;
        boundary_pair* pair;

        qgraph_edge_statistics(EdgeWeight _improvement, 
                               boundary_pair* bp, 
                               bool change) : improvement(_improvement), something_changed(change), pair(bp){
        }
};

class quotient_graph_scheduling {
        public:
                quotient_graph_scheduling();
                virtual ~quotient_graph_scheduling();

                virtual bool hasFinished() = 0;
                virtual boundary_pair & getNext() = 0;
                virtual void pushStatistics(qgraph_edge_statistics & statistic) = 0;

};


#endif /* end of include guard: QUOTIENT_GRAPH_SCHEDULING_NEFT9H3J */
