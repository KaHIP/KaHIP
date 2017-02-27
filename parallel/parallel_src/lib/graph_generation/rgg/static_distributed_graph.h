// File:   static_distributed_graph.h
// Author: Manuel Holtgrewe <holtgrewe@ira.uka.de>
//
// This header defines a static distributed graph.
//
// Types:
//
//   StaticDistributedGraph Class that represents a distributed graph, contains
//                          a local and a ghost graph.
//
// Nomenclauture
// -------------
//
// A Graph G = (V, E) has a set V of vertices (not nodes) and a set E of
// edges.
//
// Given an edge (u, v), u is called the tail, v is called the head.
//
// A "star" around a vertex u is the pair (u, { (u, v) \in E }), i.e. a
// vertex with its outgoing edges.
//
// The vertices that are owned by a process are "local vertices".  The
// vertices of which a process knows the stars but which are not owned
// by it are called "neighbour vertices".  Other vertices are called "remote
// vertices".
//
// A "splitters array" or just "splitters" is array with nondescending values.
// These splitters split a sequence of length n into k parts.  The first value
// is always 0, the last one is always n+1, there are k+1 entries.  The
// splitters array contains intervals, the i-th interval of indices in the
// sequence is [splitters[i], splitters[i+1]).
//
// Naming conventions
// ------------------
//
// edge identifier variables    e, f, g, ...
// vertex identifier variables  u, v, w, ...
//
// Vertex Ownership
// ----------------
//
// On coarsening, when the hierarchy is constructed, we give out new vertex
// numbers on every level.  Each process holds contiguously numbered vertices.
//
// Since very process knows about all its neighbour vertices, the neighbouring
// processes can send update notifications when a vertex is migrated away from
// them to a third process.  The same is true for edges.

#ifndef STATIC_DISTRIBUTED_GRAPH_H
#define STATIC_DISTRIBUTED_GRAPH_H

#include "macros_common.h"
#include "graph_types.h"
#include "util.h"

// Global Forward-Declarations.
namespace MPI { class Intracomm; }

namespace distributed_graph {

// Forward-Declarations.
class StaticLocalGraph;
class StaticLocalGraphProperties;
class StaticNeighbourVertexProperties;

// A static, distributed graph.
//
// This graph implementation does not allow migration of vertices.
class StaticDistributedGraph
{
public:
  // The default constructor initialize all members with 0.
  StaticDistributedGraph();

  // Constructor...
  ~StaticDistributedGraph();

  // Initialize the object.
  //
  // The number of processes will be taken from communicator, splitters are
  // calculated from the global values of n and m and the number of processes in
  // the communicator.
  //
  // Also compute the neighbour vertex properties.
  //
  // TODO(manuel): global_{n, m} can be inflected from splitters!
  //
  // Parameters:
  //
  //   communicator       The owning communicator.
  //   global_n           Global number of vertices.
  //   global_m           Global number of edges.
  //   local_graph        Pointer to the local graph.
  void Init(const MPI::Intracomm *communicator,
            long global_n,
            long global_m,
            const VertexId *vertex_splitters,
            const EdgeId *edge_splitters,
            const StaticLocalGraph *local_graph,
            StaticLocalGraphProperties *local_properties);

  // Returns the global number of vertices.
  const long global_n() const
  { return global_n_; }

  // Returns the global number of edges.
  const long global_m() const
  { return global_m_; }

  // Returns the local number of vertices.
  const long local_n() const
  { return local_n_; }

  // Returns the local number of edges.
  const long local_m() const
  { return local_m_; }

  // Returns the vertex splitters array.
  const VertexId *vertex_splitters() const
  {
    return vertex_splitters_;
  }

  // Returns the edge splitters array.
  const EdgeId *edge_splitters() const
  {
    return edge_splitters_;
  }

  // Returns a pointer to the local graph.
  const StaticLocalGraph *local_graph() const
  { return local_graph_; }

  // Returns a pointer to the local graph properties.
  const StaticLocalGraphProperties *local_graph_properties() const
  { return local_graph_properties_; }

  // Returns a pointer to the local graph properties, non-const version.
  StaticLocalGraphProperties *local_graph_properties()
  { return local_graph_properties_; }

  // Returns a pointer to the neighbour vertex properties object.
  const StaticNeighbourVertexProperties *neighbour_vertex_properties() const
  { return neighbour_vertex_properties_; }

  // Return the communicator in which this graph is distributed.
  const MPI::Intracomm *communicator() const
  { return &communicator_; }

  // Adds up the vertex and edge weights globally of all local graphs.
  void AddUpWeights(VertexCount *n, EdgeCount *m) const;

  // Set the local graph and the local graph properites into *g and *p and set
  // the corresponding members to NULL.
  //
  // TODO(manuel): Well, this is a hack, isn't it? Needed for LoadGraphCentrally with .daa format. Think of something better!
  void UnsafeLiberateLocalParts(StaticLocalGraph **g, StaticLocalGraphProperties **p)
  {
    *g = const_cast<StaticLocalGraph*>(local_graph_);
    local_graph_ = NULL;
    *p = const_cast<StaticLocalGraphProperties*>(local_graph_properties_);
    local_graph_properties_ = NULL;
  }

private:
  // Build neighbour_vertex_properties_.
  void InitNeighbourVertexProperties();

  // Number of vertices in the whole graph.
  long global_n_;
  // Number of edges in the whole graph.
  long global_m_;

  // Number of vertices in the local portion.
  long local_n_;
  // Number of edges in the local portion.
  long local_m_;

  // Offset of the first local vertex.
  long local_vertex_offset_;
  // Offset of the first local edge.
  long local_edge_offset_;

  // In a graph distributed on p processes, this array will have p+1 entries.
  // The vertex_splitters_[i] gives the global id of the first vertex that is
  // stored locally and vertex_splitters_[i] gives the global id of the first
  // vertex that is stored by the next process.
  const VertexId *vertex_splitters_;
  // Similar to vertex_splitters_ but for edges.
  const EdgeId *edge_splitters_;

  // The static local graph.
  const StaticLocalGraph *local_graph_;
  // The static local graph properties
  StaticLocalGraphProperties *local_graph_properties_;
  // The properties of the neighbour vertices; 
  StaticNeighbourVertexProperties *neighbour_vertex_properties_;

  // The distributed graph is distributed in this intracommunicator.
  //
  // This communicator is copied into the distributed graph on initialization.
  MPI::Intracomm communicator_;

  DISALLOW_COPY_AND_ASSIGN(StaticDistributedGraph);
};

}  // namespace distributed_graph

#endif /* ifndef STATIC_DISTRIBUTED_GRAPH_H */

