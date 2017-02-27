// File:   generate_grid.cpp
// Author: Manuel Holtgrewe <holtgrewe@ira.uka.de>

#include "generate_grid.h"  // This file's header.

#include <mpi.h>

#include "error_codes.h"
#include "graph_io.h"
#include "static_distributed_graph.h"
#include "static_local_graph.h"
#include "static_local_graph_properties.h"
#include "util.h"

namespace distributed_graph {

// Private Declarations
//

// Returns true iff the vertex at (x, y, z) has a neighbour in direction
// (dx, dy, dz) given the grid dimensions (width, height, depth);
inline
bool HasNeighbourInDirection(long x, long y, long z,
                             long dx, long dy, long dz,
                             long width, long height, long depth);


// Returns the vertex identifier for the given vertex in the given direction.
inline
VertexId NeighbourIdInDirection(long x, long y, long z,
                                long dx, long dy, long dz,
                                long width, long height, long depth);


// Implementation
//

inline
bool HasNeighbourInDirection(long x, long y, long z,
                             long dx, long dy, long dz,
                             long width, long height, long depth)
{
  ASSERT_EQ( iabs(dx) + iabs(dy) + iabs(dz), 1 );
  if (x + dx < 0 or x + dx >= width) return false;
  if (y + dy < 0 or y + dy >= height) return false;
  if (z + dz < 0 or z + dz >= depth) return false;
  return true;
}


inline
VertexId NeighbourIdInDirection(long x, long y, long z,
                                long dx, long dy, long dz,
                                long width, long height, long depth)
{
  ASSERT_TRUE( HasNeighbourInDirection(x, y, z, dx, dy, dz,
                                       width, height, depth) );
  return (x + dx) + (y + dy) * width + (z + dz) * width * height;
}


int GenerateGridGraph(long width, long height, long depth,
                          StaticDistributedGraph **graph,
                          VertexCoordinate **coordinates)
{
  ASSERT_GEQ( depth, height );
  ASSERT_GEQ( depth, width );

  EdgeId m = 6 * width * height * depth - 2 * width * height -
    2 * width * depth - 2 * height * depth;

  int process_count = MPI::COMM_WORLD.Get_size();
  int rank = MPI::COMM_WORLD.Get_rank();
  if (depth < process_count) {
    fprintf(stderr, "ERROR: Grid depth has to be at least the number"
            " of processes.\n");
    return kInvalidArguments;
  }

  // Compute vertex splitters.
  VertexId *vertex_splitters;
  ComputeSplitters(width * height * depth, process_count, &vertex_splitters);
  // Compute depth splitters, local n and and offset of the first vertex.
  VertexId *depth_splitters;
  ComputeSplitters(depth, process_count, &depth_splitters);
  long local_depth = depth_splitters[rank + 1] - depth_splitters[rank];
  VertexId local_n = width * height * local_depth;
  VertexId vertex_offset = vertex_splitters[rank];
  // Compute local m.
  EdgeCount local_m = 6 * local_n - 2 * width * local_depth -
    2 * height * local_depth;
  if (rank == 0)
    local_m -= width * height;
  if (rank == process_count - 1)
    local_m -= width * height;
  // Compute the edge splitters and edge offset.
  EdgeId *edge_splitters = new EdgeId[process_count + 1];
  edge_splitters[0] = 0;
  long n_on_0 = vertex_splitters[1] - vertex_splitters[0];
  long depth_on_0 = depth_splitters[1] - depth_splitters[0];
  edge_splitters[1] = 6 * n_on_0 - 2 * width * depth_on_0 -
    2 * height * depth_on_0 - width * height;
  for (long i = 1; i < process_count - 1; ++i) {
    long n_on_i = vertex_splitters[i + 1] - vertex_splitters[i];
    long depth_on_i = depth_splitters[i + 1] - depth_splitters[i];
    vertex_splitters[i + 1] += 6* n_on_i - 2 * width * depth_on_i -
      2 * height * depth_on_i;
  }
  edge_splitters[process_count] = m;
  EdgeId edge_offset = edge_splitters[rank];
  // Compute global n and m.
  VertexId global_n = width * height * depth;
  EdgeId global_m = 6 * width * height * depth - 2 * width * height -
    2 * width * depth - 2 * height * depth;

  // Allocate memory for adjacency datastructures.
  EdgeId *bucket_offsets = new EdgeId[local_n + 1];
  VertexId *edge_tails = new VertexId[local_m];
  VertexId *edge_heads = new VertexId[local_m];
  *coordinates = new VertexCoordinate[3 * local_n];

  // Fill adjacency array structure.
  VertexId next_u = 0;
  EdgeId next_e = 0;
  for (long z = depth_splitters[rank], zend = depth_splitters[rank + 1];
       z != zend; ++z) {
    for (long y = 0; y < height ; ++y) {
      for (long x = 0; x < width; ++x) {
        VertexId u = next_u++;
        bucket_offsets[u] = edge_offset + next_e;
        (*coordinates)[3 * u] = x;
        (*coordinates)[3 * u + 1] = y;
        (*coordinates)[3 * u + 2] = z;
        ASSERT_EQ( u, x + y * width + z * width * height );
        // Create edges for the vertex (x, y, z).
        for (long dx = -1; dx <= 1; ++dx) {
          for (long dy = -1; dy <= 1; ++dy) {
            for (long dz = -1; dz <= 1; ++dz) {
              if (iabs(dx) + iabs(dy) + iabs(dz) != 1) continue;
              if (HasNeighbourInDirection(x, y, z, dx, dy, dz, width,
                                          height, depth)) {
                EdgeId e = next_e++;
                edge_tails[e] = u;
                edge_heads[e] = NeighbourIdInDirection(x, y, z, dx, dy, dz,
                                                       width, height, depth);
              }
            }
          }
        }
      }
    }
  }
  ASSERT_EQ( next_u, local_n );
  ASSERT_EQ( next_e, local_m );
  bucket_offsets[next_u] = edge_offset + next_e;
  // Create local graph object. 
  StaticLocalGraph *local_graph = new StaticLocalGraph();
  local_graph->Init(vertex_offset, edge_offset, local_n, local_m,
                    bucket_offsets, edge_tails, edge_heads);

  // Generate the graph properties.
  StaticLocalGraphProperties *properties = new StaticLocalGraphProperties();
  properties->InitEmpty(vertex_offset, edge_offset, local_n, local_m);

  // Return the distributed graph.
  if (rank == 0)
    fprintf(stderr, "Created %ld x %ld x %ld grid with n = %lu and m = %lu\n",
            width, height, depth, global_n, global_m);
  *graph = new StaticDistributedGraph();
  (*graph)->Init(&MPI::COMM_WORLD, global_n, global_m, vertex_splitters,
                 edge_splitters, local_graph, properties);

  // Cleanup.
  delete[] depth_splitters;
  return kOk;
}

}  // namespace distributed_graph

