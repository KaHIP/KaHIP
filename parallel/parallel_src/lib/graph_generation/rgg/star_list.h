/******************************************************************************
 * star_list.h
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


#ifndef STAR_LIST_H
#define STAR_LIST_H

#include "graph_types.h"
#include "static_distributed_graph.h"
#include "static_local_graph.h"
#include "static_local_graph_properties.h"
#include "macros_assertions.h"
#include "macros_common.h"
#include "util_memory.h"


namespace distributed_graph {


// Represents a part of a static, sparse graph.  Used for sending graph
// fragments via MPI.
//
// This is a list of vertices with their outgoing edges.  Each vertex is
// assigned a temporary vertex identifier (an index in the arrays of this
// datastructure).  The StarList then consists of a mapping from temporary
// vertex id to the original vertex id and an adjacency array encoding of
// edges, the edge heads are real vertex identifiers.
//
// Additionally, the star list stores the owning process for each edge
// head.
//
// Note that this only works when VertexCount == EdgeCount == VertexId ==
// ProcessId.
class StarList
{
public:
  // The constructor initializes the StarList with null values.
  StarList() :
    n_(0), m_(0), real_vertex_ids_(NULL), edge_bucket_offsets_(NULL),
    edge_heads_(NULL)/*, edge_head_owners_(NULL)*/
  {}

  // Initialize with NULL buffer as if default-constructed.
  void Clear()
  {
    n_ = 0;
    m_ = 0;
    real_vertex_ids_ = NULL;
    edge_bucket_offsets_ = NULL;
    edge_heads_ = NULL;
  }

  // Initialize an empty star list with the given number of vertices and
  // edges on the given buffer.
  void InitEmptyOnBuffer(VertexCount n, EdgeCount m, VertexId *buffer)
  {
    n_ = n;
    m_ = m;
    buffer_ = buffer;
    buffer_[0] = n_;
    buffer_[1] = m_;

    real_vertex_ids_ = buffer + 2;
    edge_bucket_offsets_ = real_vertex_ids_ + n_;
    edge_bucket_offsets_[0] = 0;
    edge_heads_ = edge_bucket_offsets_ + n_ + 1;
    //edge_head_owners_ = edge_head_owners_ + m_;
  }

  // Initialize star list from buffer.
  //
  // Note that the star list takes NO ownership of the buffer.
  void InitFromBuffer(VertexId *buffer)
  {
    n_ = buffer[0];
    m_ = buffer[1];

    real_vertex_ids_ = buffer + 2;
    edge_bucket_offsets_ = real_vertex_ids_ + n_;
    edge_bucket_offsets_[0] = 0;
    edge_heads_ = edge_bucket_offsets_ + n_ + 1;
    //edge_head_owners_ = edge_head_owners_ + m_;
  }

  // Returns the size of a buffer as a MPI::INTEGER count for the given values
  // of n and m.
  static VertexCount BufferSizeFor(VertexCount n, EdgeCount m)
  {
    ASSERT_GEQ( n, 0 );
    ASSERT_GEQ( m, 0 );
    return 2 + n + n + 1 + m /* + m */;
  }

  VertexCount n() const
  { return n_; }

  EdgeCount m() const
  { return m_; }

  // Sets the real vertex identifier v for temporary vertex identifier i.
  void set_real_vertex_id(VertexId i, VertexId v)
  {
    ASSERT_BETWEEN( 0, i, n_ - 1 );
    ASSERT_GEQ( v, 0 );
    real_vertex_ids_[i] = v;
  }

  // Returns the real vertex identifier for temporary vertex identifier i.
  VertexId real_vertex_id(VertexId i) const
  {
    ASSERT_BETWEEN( 0, i, n_ - 1 );
    return real_vertex_ids_[i];
  }

  // Set the out degree of the vertex with temporary id i.
  //
  // Note: This function will only work correctly when called with a non-
  // decreasing value of i between the calls, skipping no value for i. Also
  // note that the d's must sum up to m.
  void set_out_degree(VertexId i, EdgeCount d)
  {
    ASSERT_BETWEEN( 0, i, n_ - 1 );
    ASSERT_BETWEEN( 0, edge_bucket_offsets_[i] + d, m_ );
    edge_bucket_offsets_[i + 1] = edge_bucket_offsets_[i] + d;
  }

  // Returns the out degree of the vertex with temporary identifier i.
  EdgeCount out_degree(VertexId i) const
  {
    ASSERT_BETWEEN( 0, i, n_ - 1 );
    return edge_bucket_offsets_[i + 1] - edge_bucket_offsets_[i];
  }


  // Set the j-th neighbour of the vertex with temporary id idx to the global
  // vertex identifier v.
  void set_neighbour(VertexId idx, EdgeCount j, VertexId v)
  {
    ASSERT_BETWEEN( 0, idx, n_ - 1 );
    ASSERT_BETWEEN( 0, j, out_degree(idx) );
    ASSERT_GEQ( v, 0 );
    edge_heads_[edge_bucket_offsets_[idx] + j] = v;
  }

  VertexId neighbour(VertexId idx, EdgeCount j) const
  {
    ASSERT_BETWEEN( 0, idx, n_ - 1 );
    ASSERT_BETWEEN( 0, j, out_degree(idx) );
    return edge_heads_[edge_bucket_offsets_[idx] + j];
  }

  /*
  // Return the owner of the j-th neighbour of idx.
  ProcessId neighbour_owner(VertexId idx, EdgeCount j) const
  {
    ASSERT_BETWEEN( 0, idx, n_ - 1 );
    ASSERT_BETWEEN( 0, j, out_degree(idx) );
    return edge_head_owners_[edge_bucket_offets_[idx] + j];
  }

  // Set the owner of the j-th neighbour of idx.
  void set_neighbour_owner(VertexId idx, EdgeCount j, ProcessId pid)
  {
    ASSERT_BETWEEN( 0, idx, n_ - 1 );
    ASSERT_BETWEEN( 0, j, out_degree(idx) );
    edge_head_owners_[edge_bucket_offets_[idx] + j] = pid;
  }
  */

private:
  // The buffer the StarList is implemented on. 
  //
  // The data layout is as follows:
  //
  //   1   * int n
  //   1   * int m
  //   n   * int owning_process_
  //   n+1 * int edge_bucket_offsets_
  //   m   * int edge_heads_
  ////   m   * int edge_head_owners_
  VertexId *buffer_;

  // Number of vertices/stars in list.
  VertexCount n_;

  // Number of edges in list.
  EdgeCount m_;

  // Mapping from temporary vertex id to real vertex id.
  VertexId *real_vertex_ids_;

  // Mapping from temporary vertex id to the owning process id.
  ProcessId *owning_process_id_;

  // Offsets into the CSR structure;
  VertexId *edge_bucket_offsets_;

  // Edge heads of the CSR structure.
  VertexId *edge_heads_;

  // Owner of the edge heads heads in the CSR structure.
  //VertexId *edge_head_owners_;

  DISALLOW_COPY_AND_ASSIGN(StarList);
};


// Similar to StarList but also sends the following vertex and edge
// attributes:
//
//   - vertex inner edges
//   - vertex weight
//   - edge weight
//
// Note that this only works when VertexCount == EdgeCount == VertexId ==
// ProcessId.
class AttributedStarList
{
public:
  // The constructor initializes the StarList with null values.
  AttributedStarList() :
    n_(0), m_(0), real_vertex_ids_(NULL), edge_bucket_offsets_(NULL),
    edge_heads_(NULL), vertex_inner_edges_(NULL), vertex_weights_(NULL),
    edge_weights_(NULL), vertex_offset_(0)
  {}

  // Initialize with NULL buffer as if default-constructed.
  void Clear()
  {
    n_ = 0;
    m_ = 0;
    real_vertex_ids_ = NULL;
    edge_bucket_offsets_ = NULL;
    edge_heads_ = NULL;
    vertex_inner_edges_ = NULL;
    vertex_weights_ = NULL;
    edge_weights_ = NULL;
    vertex_offset_ = 0;
  }

  // Initialize an empty star list with the given number of vertices and
  // edges on the given buffer.
  void InitEmptyOnBuffer(VertexCount n, EdgeCount m, VertexId *buffer)
  {
    n_ = n;
    m_ = m;
    buffer_ = buffer;
    buffer_[0] = n_;
    buffer_[1] = m_;

    real_vertex_ids_ = buffer_ + 2;
    edge_bucket_offsets_ = real_vertex_ids_ + n_;
    edge_bucket_offsets_[0] = 0;
    edge_heads_ = edge_bucket_offsets_ + n_ + 1;
    vertex_inner_edges_ = edge_heads_ + m_;
    vertex_weights_ = vertex_inner_edges_ + n_;
    edge_weights_ = vertex_weights_ + n_;

    vertex_offset_ = 0;
  }

  // Initialize star list from buffer.
  //
  // Note that the star list takes NO ownership of the buffer.
  void InitFromBuffer(VertexId *buffer)
  {
    buffer_ = buffer;
    n_ = buffer_[0];
    m_ = buffer_[1];

    real_vertex_ids_ = buffer_ + 2;
    edge_bucket_offsets_ = real_vertex_ids_ + n_;
    ASSERT_EQ( edge_bucket_offsets_[0], 0 );
    ASSERT_EQ( edge_bucket_offsets_[n_], m_ );
    edge_heads_ = edge_bucket_offsets_ + n_ + 1;
    vertex_inner_edges_ = edge_heads_ + m_;
    vertex_weights_ = vertex_inner_edges_ + n_;
    edge_weights_ = vertex_weights_ + n_;
  }

  // Returns the size of a buffer as a MPI::INTEGER count for the given values
  // of n and m.
  static VertexCount BufferSizeFor(VertexCount n, EdgeCount m)
  {
    ASSERT_GEQ( n, 0 );
    ASSERT_GEQ( m, 0 );
    return 2 + n + n + 1 + m + n + n + m;
  }

  // Copy in the star around v (with attributes).
  void CopyStarIn(const StaticDistributedGraph *graph, const VertexId v)
  {
    // Get shortcuts.
    const StaticLocalGraphProperties *properties =
      graph->local_graph_properties();
    const StaticLocalGraph *local_graph = graph->local_graph();
    // Copy in the vertices' properties.
    set_real_vertex_id(vertex_offset_, v);
    set_vertex_weight(vertex_offset_, properties->vertex_weight(v));
    set_vertex_inner_edges(vertex_offset_, properties->internal_edges(v));
    set_out_degree(vertex_offset_, local_graph->out_degree(v));
    // Copy in the edges and their weight.
    for (int i = 0, iend = local_graph->out_degree(v); i < iend; ++i) {
      set_neighbour(vertex_offset_, i, local_graph->neighbour(v, i));
      EdgeId e = local_graph->out_edge(v, i);
      //printf("copy in weight: %d\n", properties->edge_weight(e)); // XXX
      set_edge_weight(vertex_offset_, i, properties->edge_weight(e));
    }
    vertex_offset_ += 1;
  }

  VertexCount n() const
  { return n_; }

  EdgeCount m() const
  { return m_; }

  // Sets the real vertex identifier v for temporary vertex identifier i.
  void set_real_vertex_id(VertexId i, VertexId v)
  {
    ASSERT_BETWEEN( 0, i, n_ - 1 );
    ASSERT_GEQ( v, 0 );
    real_vertex_ids_[i] = v;
  }

  // Returns the real vertex identifier for temporary vertex identifier i.
  VertexId real_vertex_id(VertexId i) const
  {
    ASSERT_BETWEEN( 0, i, n_ - 1 );
    return real_vertex_ids_[i];
  }

  // Set the out degree of the vertex with temporary id i.
  //
  // Note: This function will only work correctly when called with a non-
  // decreasing value of i between the calls, skipping no value for i. Also
  // note that the d's must sum up to m.
  void set_out_degree(VertexId i, EdgeCount d)
  {
    ASSERT_BETWEEN( 0, i, n_ - 1 );
    ASSERT_BETWEEN( 0, edge_bucket_offsets_[i] + d, m_ );
    edge_bucket_offsets_[i + 1] = edge_bucket_offsets_[i] + d;
  }

  // Returns the out degree of the vertex with temporary identifier i.
  EdgeCount out_degree(VertexId i) const
  {
    ASSERT_BETWEEN( 0, i, n_ - 1 );
    return edge_bucket_offsets_[i + 1] - edge_bucket_offsets_[i];
  }


  // Set the j-th neighbour of the vertex with temporary id idx to the global
  // vertex identifier v.
  void set_neighbour(VertexId idx, EdgeCount j, VertexId v)
  {
    ASSERT_BETWEEN( 0, idx, n_ - 1 );
    ASSERT_BETWEEN( 0, j, out_degree(idx) );
    ASSERT_GEQ( v, 0 );
    edge_heads_[edge_bucket_offsets_[idx] + j] = v;
  }

  VertexId neighbour(VertexId idx, EdgeCount j) const
  {
    ASSERT_BETWEEN( 0, idx, n_ - 1 );
    ASSERT_BETWEEN( 0, j, out_degree(idx) );
    return edge_heads_[edge_bucket_offsets_[idx] + j];
  }

  // Set the number of inner edges for the given vertex.
  void set_vertex_inner_edges(VertexId idx, EdgeCount c)
  {
    ASSERT_BETWEEN( 0, idx, n_ - 1 );
    ASSERT_GEQ( c, 0 );
    vertex_inner_edges_[idx] = c;
  }

  // Returns the number of inner edges for the given vertex.
  VertexId vertex_inner_edges(VertexId idx) const
  {
    ASSERT_BETWEEN( 0, idx, n_ - 1 );
    return vertex_inner_edges_[idx];
  }

  // Set the weight of the given vertex.
  void set_vertex_weight(VertexId idx, VertexWeight w)
  {
    ASSERT_BETWEEN( 0, idx, n_ - 1 );
    ASSERT_GT( w, 0 );
    vertex_weights_[idx] = w;
  }

  // Returns the weight of the given vertex.
  VertexId vertex_weight(VertexId idx) const
  {
    ASSERT_BETWEEN( 0, idx, n_ - 1 );
    return vertex_weights_[idx];
  }

  // Set the weight of the edge to the j-th neighbour of the vertex with
  // temporary id idx to the global vertex identifier v.
  void set_edge_weight(VertexId idx, EdgeCount j, EdgeWeight w)
  {
    ASSERT_BETWEEN( 0, idx, n_ - 1 );
    ASSERT_BETWEEN( 0, j, out_degree(idx) );
    ASSERT_GT( w, 0 );
    //printf("edge_weights_[%d] = %d\n", edge_bucket_offsets_[idx] + j, w); // XXX
    edge_weights_[edge_bucket_offsets_[idx] + j] = w;
  }

  VertexId edge_weight(VertexId idx, EdgeCount j) const
  {
    ASSERT_BETWEEN( 0, idx, n_ - 1 );
    ASSERT_BETWEEN( 0, j, out_degree(idx) );
    return edge_weights_[edge_bucket_offsets_[idx] + j];
  }

private:
  // The buffer the StarList is implemented on. 
  //
  // The data layout is as follows:
  //
  //   1   * int n
  //   1   * int m
  //   n   * int owning_process_
  //   n+1 * int edge_bucket_offsets_
  //   m   * int edge_heads_
  //   m   * int edge_head_owners_
  //   n   * int vertex inner edges
  //   n   * int vertex weights_
  //   m   * int edge_weights_
  VertexId *buffer_;

  // Number of vertices/stars in list.
  VertexCount n_;

  // Number of edges in list.
  EdgeCount m_;

  // Mapping from temporary vertex id to real vertex id.
  VertexId *real_vertex_ids_;

  // Mapping from temporary vertex id to the owning process id.
  ProcessId *owning_process_id_;

  // Offsets into the CSR structure;
  VertexId *edge_bucket_offsets_;

  // Edge heads of the CSR structure.
  VertexId *edge_heads_;

  // Number of inner edges of the vertices.
  EdgeCount *vertex_inner_edges_;

  // Vertex weights.
  VertexWeight *vertex_weights_;

  // The edge weights.
  EdgeWeight *edge_weights_;

  // The number of stars copied in via CopyStarIn().
  VertexCount vertex_offset_;

  DISALLOW_COPY_AND_ASSIGN(AttributedStarList);
};


// An UncoarsenStarList that contains the following information about the
// stars:
//
//   - real vertex ids
//   - compressed vertices
//
// It is used to uncoarsen coarser vertices into finer vertices again.
// The naming is somewhat bad, it is no real "star list".
class UncoarsenStarList
{
public:
  // The constructor initializes the StarList with null values.
  UncoarsenStarList() :
    n_(0), real_vertex_ids_(NULL),
    compressed_vertices_count_(0), compressed_vertices_offsets_(NULL),
    compressed_vertices_(NULL)
  {}

  // Initialize an empty star list with the given number of vertices and
  // Initialize with NULL buffer as if default-constructed.
  void Clear()
  {
    n_ = 0;
    real_vertex_ids_ = NULL;
    compressed_vertices_count_ = 0;
    compressed_vertices_offsets_ = NULL;
    compressed_vertices_ = NULL;
  }

  // edges on the given buffer.
  void InitEmptyOnBuffer(VertexCount n,
                         VertexCount compressed_vertices_count,
                         VertexId *buffer)
  {
    n_ = n;
    compressed_vertices_count_ = compressed_vertices_count;
    buffer_ = buffer;
    buffer_[0] = n_;
    buffer_[1] = compressed_vertices_count_;

    real_vertex_ids_ = buffer + 2;

    compressed_vertices_offsets_ = real_vertex_ids_ + n_;
    compressed_vertices_offsets_[0] = 0;
    compressed_vertices_ = compressed_vertices_offsets_ + n_ + 1;
  }

  // Initialize star list from buffer.
  //
  // Note that the star list takes NO ownership of the buffer.
  void InitFromBuffer(VertexId *buffer)
  {
    n_ = buffer[0];
    compressed_vertices_count_ = buffer[1];

    real_vertex_ids_ = buffer + 2;
    compressed_vertices_offsets_ = real_vertex_ids_ + n_;
    ASSERT_EQ( compressed_vertices_offsets_[0], 0 );
    compressed_vertices_ = compressed_vertices_offsets_ + n_ + 1;
  }

  // Returns the size of a buffer as a MPI::INTEGER count for the given values
  // of n and m.  Compressed vertex count is c.
  static VertexCount BufferSizeFor(VertexCount n, VertexCount c)
  {
    ASSERT_GEQ( n, 0 );
    ASSERT_GEQ( c, 0 );
    return 2 + n + n + 1 + c;
  }

  VertexCount n() const
  { return n_; }

  VertexCount compressed_vertices_count() const
  { return compressed_vertices_count_; }

  // Sets the real vertex identifier v for temporary vertex identifier i.
  void set_real_vertex_id(VertexId i, VertexId v)
  {
    ASSERT_BETWEEN( 0, i, n_ - 1 );
    ASSERT_GEQ( v, 0 );
    real_vertex_ids_[i] = v;
  }

  // Returns the real vertex identifier for temporary vertex identifier i.
  VertexId real_vertex_id(VertexId i) const
  {
    ASSERT_BETWEEN( 0, i, n_ - 1 );
    return real_vertex_ids_[i];
  }

  // Set the number of compressed vertices for the vertex with index i to c.
  //
  // Note that you have to set the count for i before you can set the count
  // for vertex i+1 and before you can set any compressed vertex for i.
  void set_compressed_vertex_count(VertexId i, VertexCount c)
  {
    ASSERT_BETWEEN( 0, i, n_ - 1 );
    ASSERT_BETWEEN( 0, compressed_vertices_offsets_[i] + c, compressed_vertices_count_ );
    compressed_vertices_offsets_[i + 1] = compressed_vertices_offsets_[i] + c;
  }

  // Returns the number of compressed vertices for the vertex with id i.
  VertexCount compressed_vertices_count(VertexId i) const
  {
    ASSERT_BETWEEN( 0, i, n_ - 1 );
    return compressed_vertices_offsets_[i + 1] - compressed_vertices_offsets_[i];
  }

  // Sets the j-th compressed vertex for the vertex with the index i to v.
  void set_compressed_vertex(VertexId i, VertexId j, VertexId v)
  {
    ASSERT_BETWEEN( 0, i, n_ - 1 );
    ASSERT_BETWEEN( 0, j, compressed_vertices_count(i) );
    ASSERT_GEQ( v, 0 );
    compressed_vertices_[compressed_vertices_offsets_[i] + j] = v;
  }

  // Return the identifier of the j-th compressed vertex for the vetex with id i.
  VertexId compressed_vertex(VertexId i, VertexId j) const
  {
    ASSERT_BETWEEN( 0, i, n_ - 1 );
    ASSERT_BETWEEN( 0, j, compressed_vertices_count(i) );
    return compressed_vertices_[compressed_vertices_offsets_[i] + j];
  }

private:
  // The buffer the StarList is implemented on. 
  //
  // The data layout is as follows:
  //
  //   1   * int n
  //   1   * int compressed_vertex_count_ = c
  //   n   * int real_vertex_id_
  //   n+1 * int next_vertices_offsets
  //   c   * int next_vertices
  VertexId *buffer_;

  // Number of vertices/stars in list.
  VertexCount n_;

  // Mapping from temporary vertex id to real vertex id.
  VertexId *real_vertex_ids_;

  // Number of compressed vertices.
  VertexCount compressed_vertices_count_;

  // Offset-array pointer for CSR-style list of compressed vertices.
  VertexId *compressed_vertices_offsets_;
 
  // Start pointer for CSR-style list of compressed vertices.
  VertexId *compressed_vertices_;

  DISALLOW_COPY_AND_ASSIGN(UncoarsenStarList);
};


// An RefinementStarList that contains the following information about the
// stars:
//
//   - real vertex ids
//   - vertex weights
//   - edge weights
//
// It is used to exchange the stars with the information necessary for
// Christian's refinement code.
class RefinementStarList
{
public:
  // The constructor initializes the StarList with null values.
  RefinementStarList() :
    n_(0), m_(0), real_vertex_ids_(NULL), vertex_weights_(NULL),
    edge_bucket_offsets_(NULL), edge_heads_(NULL), edge_weights_(NULL)
  {}

  // Initialize with NULL buffer as if default-constructed.
  void Clear()
  {
    n_ = 0;
    m_ = 0;
    real_vertex_ids_ = NULL;
    vertex_weights_ = NULL;
    edge_bucket_offsets_ = NULL;
    edge_heads_ = NULL;
    edge_weights_ = NULL;
  }

  // Initialize an empty star list with the given number of vertices and
  // edges on the given buffer.
  void InitEmptyOnBuffer(VertexCount n, EdgeCount m,
                         VertexId *buffer)
  {
    n_ = n;
    m_ = m;
    buffer_ = buffer;
    buffer_[0] = n_;
    buffer_[1] = m_;

    real_vertex_ids_ = buffer + 2;
    vertex_weights_ = real_vertex_ids_ + n_;
    edge_bucket_offsets_ = vertex_weights_ + n_;
    edge_bucket_offsets_[0] = 0;
    edge_heads_ = edge_bucket_offsets_ + n_ + 1;
    edge_weights_ = edge_heads_ + m_;
  }

  // Initialize star list from buffer.
  //
  // Note that the star list takes NO ownership of the buffer.
  void InitFromBuffer(VertexId *buffer)
  {
    n_ = buffer[0];
    m_ = buffer[1];

    real_vertex_ids_ = buffer + 2;
    vertex_weights_ = real_vertex_ids_ + n_;
    edge_bucket_offsets_ = vertex_weights_ + n_;
    ASSERT_EQ( edge_bucket_offsets_[0], 0 );
    edge_heads_ = edge_bucket_offsets_ + n_ + 1;
    edge_weights_ = edge_heads_ + m_;
  }

  // Returns the size of a buffer as a MPI::INTEGER count for the given values
  // of n and m.  Compressed vertex count is c.
  static VertexCount BufferSizeFor(VertexCount n, EdgeCount m)
  {
    ASSERT_GEQ( n, 0 );
    ASSERT_GEQ( m, 0 );
    return 2 + n + n + n + 1 + m + m;
  }

  VertexCount n() const
  { return n_; }

  EdgeCount m() const
  { return m_; }

  // Sets the real vertex identifier v for temporary vertex identifier i.
  void set_real_vertex_id(VertexId i, VertexId v)
  {
    ASSERT_BETWEEN( 0, i, n_ - 1 );
    ASSERT_GEQ( v, 0 );
    real_vertex_ids_[i] = v;
  }

  // Returns the real vertex identifier for temporary vertex identifier i.
  VertexId real_vertex_id(VertexId i) const
  {
    ASSERT_BETWEEN( 0, i, n_ - 1 );
    return real_vertex_ids_[i];
  }

  // Sets the vertex weight of the vertex with id i.
  void set_vertex_weight(VertexId i, VertexWeight w)
  {
    ASSERT_BETWEEN( 0, i, n_ - 1 );
    ASSERT_GT( w, 0 );
    vertex_weights_[i] = w;
  }

  // Returns the real vertex identifier for temporary vertex identifier i.
  VertexWeight vertex_weight(VertexId i) const
  {
    ASSERT_BETWEEN( 0, i, n_ - 1 );
    return vertex_weights_[i];
  }

  // Set the out degree of the vertex with temporary id i.
  //
  // Note: This function will only work correctly when called with a non-
  // decreasing value of i between the calls, skipping no value for i. Also
  // note that the d's must sum up to m.
  void set_out_degree(VertexId i, EdgeCount d)
  {
    ASSERT_BETWEEN( 0, i, n_ - 1 );
    ASSERT_BETWEEN( 0, edge_bucket_offsets_[i] + d, m_ );
    edge_bucket_offsets_[i + 1] = edge_bucket_offsets_[i] + d;
  }

  // Returns the out degree of the vertex with temporary identifier i.
  EdgeCount out_degree(VertexId i) const
  {
    ASSERT_BETWEEN( 0, i, n_ - 1 );
    return edge_bucket_offsets_[i + 1] - edge_bucket_offsets_[i];
  }

  // Set the j-th neighbour of the vertex with temporary id idx to the global
  // vertex identifier v.
  void set_neighbour(VertexId idx, EdgeCount j, VertexId v)
  {
    ASSERT_BETWEEN( 0, idx, n_ - 1 );
    ASSERT_BETWEEN( 0, j, out_degree(idx) );
    ASSERT_GEQ( v, 0 );
    edge_heads_[edge_bucket_offsets_[idx] + j] = v;
  }

  VertexId neighbour(VertexId idx, EdgeCount j) const
  {
    ASSERT_BETWEEN( 0, idx, n_ - 1 );
    ASSERT_BETWEEN( 0, j, out_degree(idx) );
    return edge_heads_[edge_bucket_offsets_[idx] + j];
  }

  // Set the weight of the edge to the j-th neighbour of the vertex with
  // temporary id idx to the weight w.
  void set_edge_weight(VertexId idx, EdgeCount j, VertexWeight w)
  {
    ASSERT_BETWEEN( 0, idx, n_ - 1 );
    ASSERT_BETWEEN( 0, j, out_degree(idx) );
    ASSERT_GT( w, 0 );
    edge_weights_[edge_bucket_offsets_[idx] + j] = w;
  }

  // Return weight of edge with vertex index idx to the j-th neighbour.
  VertexId edge_weight(VertexId idx, EdgeCount j) const
  {
    ASSERT_BETWEEN( 0, idx, n_ - 1 );
    ASSERT_BETWEEN( 0, j, out_degree(idx) );
    return edge_weights_[edge_bucket_offsets_[idx] + j];
  }

private:
  // The buffer the StarList is implemented on. 
  //
  // The data layout is as follows:
  //
  //   1   * int n
  //   1   * int m
  //   n   * int real_vertex_id_
  //   n   * int vertex_weights_
  //   n+1 * int edge_bucket_offsets_
  //   m   * int edge_heads_
  //   m   * int edge_weights_
  VertexId *buffer_;

  // Number of vertices/stars in list.
  VertexCount n_;

  // Number of edges in list.
  EdgeCount m_;

  // Mapping from temporary vertex id to real vertex id.
  VertexId *real_vertex_ids_;

  // Weights of the vertices.
  VertexWeight *vertex_weights_;

  // Offsets into the CSR structure;
  VertexId *edge_bucket_offsets_;

  // Edge heads of the CSR structure.
  VertexId *edge_heads_;

  // Weights of the edges.
  EdgeWeight *edge_weights_;

  DISALLOW_COPY_AND_ASSIGN(RefinementStarList);
};


}  // namespace distributed_graph


#endif /* ifndef STAR_LIST_H */

