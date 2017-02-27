/******************************************************************************
 * generate_grid.h
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

#ifndef GENERATE_GRID_H
#define GENERATE_GRID_H

#include "graph_types.h"

namespace distributed_graph {

// Forward-Declarations.
class StaticLocalGraph;
class StaticLocalGraphProperties;
class StaticDistributedGraph;


// Interface.

// Generates a graph representing a grid with coordinates distributedly.
//
// The graph will be local to MPI::COMM_WORLD.
//
// Note that the distribution of the graph will not be perfect since each
// process generates a slice, cuts go through the depth, along width
// and height.
int GenerateGridGraph(long width, long height, long depth,
                          StaticDistributedGraph **graph,
                          VertexCoordinate **coordinates);

}  // namespace distributed_graph

#endif // ifndef GENERATE_GRID_H

