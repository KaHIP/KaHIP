// File:   graph_io.h
// Author: Manuel Holtgrewe <holtgrewe@ira.uka.de>
//
// Routines and supporting code for loading and saving graphs.
//
// Functions:
//
//   LoadGraphCentrallyAndDistribute()
//                           Load a graph centrally and then distribute it.
//   LoadGraphDistributed()  Load a graph in a completely distributed manner.
//
//
// Binary Distributed File Format
// ------------------------------
//
// The file extension is ".daa".
//
// In order to read in a graph in parallel using MPI-2's parallel I/O routines,
// we have to use a binary format.  It is structured as follows:
//
//   Datatype   Count   Description
//
//   MPI::INT *     3   Magic number, 'DAA' for Distributed Adjacency Array.
//   MPI::INT *     1   Format version, current version is 1.
//   MPI::INT *     1   Specifies whether or not to use weights, see below.
//   MPI::INT *     1   Vertex count n.
//   MPI::INT *     1   Edge count m.
//   MPI::INT *  (n+1)  Offsets in edge head array.
//   MPI::INT *     m   Edge heads.
//   MPI::INT *     n   Vertex weights, if enabled in weight flag.
//   MPI::INT *     m   Edge weights, if enabled in weight flag.
//
// The weight flag can have the following values:
//
//   00   No weights.
//   01   Edge weights only.
//   10   Vertex weights only.
//   11   Both vertex and edge weights.

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

