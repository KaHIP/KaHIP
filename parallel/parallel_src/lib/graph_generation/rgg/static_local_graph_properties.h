/******************************************************************************
 * static_local_graph_properties.h
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

#ifndef STATIC_LOCAL_GRAPH_PROPERTIES_H
#define STATIC_LOCAL_GRAPH_PROPERTIES_H

#include "static_local_graph_properties.h"  // This file's header.

#include <stdint.h>

#include "graph_types.h"
#include "macros_assertions.h"
#include "macros_common.h"


namespace distributed_graph {

// Forward-Declarations.
class StaticLocalGraph;


// Properties for a StaticLocalGraph.
class StaticLocalGraphProperties
{
public:

  // The default constructor initializes all members with 0.
  StaticLocalGraphProperties();

  // The default destructor frees the members again.
  ~StaticLocalGraphProperties();

  // Parameters:
  //
  //   vertex_offset  Offset of the first vertex identifier.
  //   edge_offset    Offset of the first edge identifier.
  //   n              Number of local vertices.
  //   m              Number of local edges.
  //   edges          An array with m edges.
  //
  // TODO(manuel): Allow for disabling certain weights to save memory.
  //               For example, internal edges could always return 0, default
  //               vertex and edge weight could be 0.
  void InitEmpty(const VertexId vertex_offset,
                 const EdgeId edge_offset,
                 const long n,
                 const long m);


  // Initialize with mapping from identifiers to weights/edge counts.
  //
  // Note that the graph object takes ownership of the arrays and will
  // deallocate them when destructed.
  //
  // Parameters:
  //
  //   vertex_offset   Offset of the first vertex identifier.
  //   edge_offset     Offset of the first edge identifier.
  //   n               Number of local vertices.
  //   m               Number of local edges.
  //   edges           An array with m edges.
  //   vertex_weights  Array with vertex weights.
  //   internal_edges  Array with internal edges counts.
  //   edge_weights    Array with edge weights.
  //
  // TODO(manuel): Allow for disabling certain weights to save memory.
  //               For example, internal edges could always return 0, default
  //               vertex and edge weight could be 0.
  void InitWithValues(const VertexId vertex_offset,
                      const EdgeId edge_offset,
                      const long n,
                      const long m,
                      VertexWeight *vertex_weights,
                      EdgeCount *internal_edges,
                      RggEdgeWeight *edge_weights);

  // Set the coordinates.
  //
  // Note that this function takes ownership of the coordinates.
  //
  // If there already are coordinates, they are freed.
  void InitCoordinates(VertexCoordinate *coordinates);

  // Access for vertex weights.
  const VertexWeight &vertex_weight(const VertexId &v) const
  {
    ASSERT_BETWEEN( vertex_offset_, v, vertex_offset_ + n_ - 1);
    VertexId vl = v - vertex_offset_;
    return vertex_weights_[vl];
  }

  // Mutator for vertex weights.
  void set_vertex_weight(const VertexId &v, const VertexWeight &w)
  {
    ASSERT_BETWEEN( vertex_offset_, v, vertex_offset_ + n_ - 1);
    VertexId vl = v - vertex_offset_;
    vertex_weights_[vl] = w;
  }

  // Access for internal edge count.
  const EdgeCount &internal_edges(const VertexId &v) const
  {
    ASSERT_BETWEEN( vertex_offset_, v, vertex_offset_ + n_ - 1);
    VertexId vl = v - vertex_offset_;
    return internal_edges_[vl];
  }

  // Mutator for internal edge counts.
  void set_internal_edges(const VertexId &v, const EdgeCount &c)
  {
    ASSERT_BETWEEN( vertex_offset_, v, vertex_offset_ + n_ - 1);
    VertexId vl = v - vertex_offset_;
    internal_edges_[vl] = c;
  }

  // Access for edge weights.
  const VertexWeight &edge_weight(const EdgeId &e) const
  {
    ASSERT_BETWEEN( edge_offset_, e, edge_offset_ + m_ -1);
    EdgeId local_e = e - edge_offset_;
    return edge_weights_[local_e];
  }

  // Mutator for edge weights.
  void set_edge_weight(const EdgeId &e, const RggEdgeWeight &w)
  {
    ASSERT_BETWEEN( edge_offset_, e, edge_offset_ + m_ - 1 );
    EdgeId el = e - edge_offset_;
    edge_weights_[el] = w;
  }

  // Access for the coordinates of vertex v.
  void coordinates(const VertexId &v, VertexCoordinate *x, VertexCoordinate *y,
                   VertexCoordinate *z) const
  {
    ASSERT_TRUE( has_coordinates() );
    ASSERT_BETWEEN( 0, v, n_ - 1 );
    *x = coordinates_[v * 3];
    *y = coordinates_[v * 3 + 1];
    *z = coordinates_[v * 3 + 2];
  }

  const VertexCoordinate *vertex_coordinates() const
  { return coordinates_; }

  // Returns true iff coordinates have been set.
  bool has_coordinates() const
  { return coordinates_ != NULL; }

  // Returns pointer to edge weights.
  const RggEdgeWeight *edge_weights() const
  { return edge_weights_; }

private:
  // Number of vertices.
  long n_;

  // Number of edges.
  long m_;

  // The vertex offset of the local graph.
  long vertex_offset_;

  // The edge offset of the local graph.
  long edge_offset_;

  // The vertex weights.
  VertexWeight *vertex_weights_;

  // The internal edge counts.
  EdgeCount *internal_edges_;

  // The edge weights.
  RggEdgeWeight *edge_weights_;

  // The vertex coordinates.  Initially, it is NULL and has to be set with
  // InitVertexCoordinates().
  VertexCoordinate *coordinates_;

  // For simplicity and performance, we allow friend-access for I/O routines.
  friend int SaveGraphCentrally(
      const GraphFileFormat format,
      const char *filename,
      const StaticLocalGraph *graph,
      const StaticLocalGraphProperties *graph_properties);

  friend class UnsafeStaticLocalGraphPropertiesAccess;

  DISALLOW_COPY_AND_ASSIGN(StaticLocalGraphProperties);
};


// This class allows unsafe access to the properties of
//
// Note that it is your responsibility not to change any of the data or to
// live with undefined behaviour.
class UnsafeStaticLocalGraphPropertiesAccess
{
public:
  // Initialize the object with a pointer to the given graph.
  UnsafeStaticLocalGraphPropertiesAccess(const StaticLocalGraphProperties *p) :
    properties_(p)
  {}

  // Return a non-const version of the vertex weights.
  VertexWeight *UNSAFE_vertex_weights()
  { return const_cast<VertexWeight*>(properties_->vertex_weights_); }

  // Return a non-const version of the vertex weights.
  EdgeCount *UNSAFE_internal_edges()
  { return const_cast<EdgeCount*>(properties_->internal_edges_); }

  // Return a non-const version of the vertex weights.
  RggEdgeWeight *UNSAFE_edge_weights()
  { return const_cast<RggEdgeWeight*>(properties_->edge_weights_); }

private:
  // A pointer to the static local graph.
  const StaticLocalGraphProperties *properties_;

  DISALLOW_COPY_AND_ASSIGN(UnsafeStaticLocalGraphPropertiesAccess);
};


}  // namespace distributed_graph


#endif /* ifndef STATIC_LOCAL_GRAPH_PROPERTIES_H */

