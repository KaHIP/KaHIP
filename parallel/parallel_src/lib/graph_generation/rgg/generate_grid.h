// File:   generate_grid.h
// Author: Manuel Holtgrewe <holtgrewe@ira.uka.de>
//
// Code to generate grid graphs.

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

