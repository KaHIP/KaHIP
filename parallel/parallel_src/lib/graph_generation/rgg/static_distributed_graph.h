/******************************************************************************
 * static_distributed_graph.h
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

