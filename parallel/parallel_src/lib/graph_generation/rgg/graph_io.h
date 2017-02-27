/******************************************************************************
 * graph_io.h
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

#ifndef GRAPH_IO_H
#define GRAPH_IO_H

#include <stdint.h>

#include <mpi.h>

#include "graph_types.h"

namespace distributed_graph {

// Forward-Declarations.
//
class StaticDistributedGraph;
class StaticLocalGraph;
class StaticLocalGraphProperties;

int SaveGraphDistributedlyMetis(const char *filename,
                               const StaticDistributedGraph *graph);


// Store the coordinates in binary format to filename distributedly.
int StoreCoordinatesDistributedlyMetis(const char *filename,
                                  const StaticDistributedGraph *graph,
                                  const VertexCoordinate *coordinates);

}  // namespace distributed_graph

#endif /* ifndef GRAPH_IO_H */

