/******************************************************************************
 * parallel_projection.h
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

#ifndef PARALLEL_PROJECTION_HBRCPQ0P
#define PARALLEL_PROJECTION_HBRCPQ0P

#include "data_structure/parallel_graph_access.h"

class parallel_projection {
public:
        parallel_projection();
        virtual ~parallel_projection();

        void parallel_project( MPI_Comm communicator, parallel_graph_access & finer, parallel_graph_access & coarser );

        //initial assignment after initial partitioning
        void initial_assignment( parallel_graph_access & G, complete_graph_access & Q);
private:
        std::vector< std::vector< NodeID > > m_messages;
};




#endif /* end of include guard: PARALLEL_PROJECTION_HBRCPQ0P */
