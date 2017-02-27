/******************************************************************************
 * static_local_graph.h
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

#ifndef STATIC_LOCAL_GRAPH_H
#define STATIC_LOCAL_GRAPH_H


#include "macros_assertions.h"
#include "macros_common.h"
#include "graph_types.h"
#include "util.h"


namespace distributed_graph {

// Forward-Declarations.
class StaticLocalGraph;  // For friend below.
class StaticLocalGraphProperties;  // For friend below.


// The local part of a distributed graph.
//
// It uses a simple adjacency array representation.  The vertex identifiers
// must be contiguous.
// TODO(manuel): At the moment, we store edge and vertex identifiers with global ids.  Maybe it would be worthwile storing only local identifiers.
class StaticLocalGraph
{
public:
  // The default constructor initializes all members of the class to 0.
  StaticLocalGraph();

  // The default destructor.
  ~StaticLocalGraph();

  // (Re-)initialize the graph with the given edge list.
  //
  // The edge list must be sorted lexicographically in ascending order:  Order
  // by tail and break ties by comparing the head of the edge.
  //
  // Runs in time O(m).
  //
  // Parameters:
  //
  //   vertex_offset  Offset of the first vertex identifier.
  //   edge_offset    Offset of the first edge identifier.
  //   n              Number of local vertices.
  //   m              Number of local edges.
  //   edges          An array with m edges.
  void Init(const VertexId vertex_offset,
            const EdgeId edge_offset,
            const VertexId n,
            const VertexId m,
            const Edge *edges);

  // (Re-)initialize the graph with the given adjacency array structure.
  //
  // Note that the graph takes ownership of the edge_offsets, edge_tails
  // and edge_heads arrays.
  //
  // Runs in time O(m).
  //
  // Parameters:
  //
  //   vertex_offset  Offset of the first vertex identifier.
  //   edge_offset    Offset of the first edge identifier.
  //   n              Number of local vertices.
  //   m              Number of local edges.
  //   bucket_begin_offsets
  //                  Array with n+1 offsets into the edge array.
  //   edge_tails     Tails of the m edges.
  //   edge_heads     Heads of the m edges.
  void Init(const VertexId vertex_offset,
            const EdgeId edge_offset,
            const VertexId n,
            const EdgeId m,
            EdgeId *bucket_begin_offsets,
            VertexId *edge_tails,
            VertexId *edge_heads);

  // Dump the graph to out.
  void Dump(FILE *out) const;

  // Dump a histogram of vertex degrees to the given FILE*.
  void DumpDegreeHistogram(FILE *out) const;

  // Returns the number of vertices.
  VertexId n() const
  { return n_; }

  // Returns the number of edges.
  VertexId m() const
  { return m_; }

  // Returns the out degree of the given vertex.
  //
  // Parameters:
  //
  //   v  Global vertex id.
  VertexId out_degree(const VertexId &v) const
  {
    VertexId local_v = v - vertex_offset_;
    ASSERT_BETWEEN( 0, local_v, n_ - 1 );
/* TODO(manuel): The following is not true when used in Manne-Bisseling matching.
#ifndef NDEBUG
    if (bucket_begin_offsets_[local_v] != bucket_begin_offsets_[local_v + 1]) {
      ASSERT_EQ( edge_tails_[bucket_begin_offsets_[local_v] - edge_offset_], v );
      if (local_v + 1 < n_)
        ASSERT_GT( edge_tails_[bucket_begin_offsets_[local_v + 1] - edge_offset_], v );
    }
#endif
*/
    VertexId res = bucket_begin_offsets_[local_v + 1] - bucket_begin_offsets_[local_v];
    ASSERT_GEQ( res, 0 );
    return res;
  }

  // Returns the global identifier of the ith neighbour of the given vertex.
  //
  // Parameters:
  //
  //   v  Global vertex id.
  //   i  Index of the neighbour to retrieve.
  VertexId neighbour(const VertexId &v, EdgeCount i) const
  {
    VertexId local_v = v - vertex_offset_;
    ASSERT_BETWEEN( 0, local_v, n_ - 1 );
    ASSERT_LT( bucket_begin_offsets_[local_v] + i,
                 bucket_begin_offsets_[local_v + 1] );
    EdgeId local_e = bucket_begin_offsets_[local_v] - edge_offset_ + i;
    ASSERT_LT( local_e, m_ );
    return edge_heads_[local_e];
  }

  // Returns the local id of the edge to the ith neighbour.
  //
  // Parameters:
  //
  //   v  Global vertex id.
  //   i  Index of the edge to retrieve.
  EdgeId out_edge(const VertexId &v, EdgeCount i) const
  {
    //std::cout <<  v       << std::endl;
    VertexId local_v = v - vertex_offset_;
    ASSERT_BETWEEN( 0, local_v, n_ - 1 );
    ASSERT_LT( i, bucket_begin_offsets_[local_v + 1] );
    EdgeId local_e = bucket_begin_offsets_[local_v] - edge_offset_ + i;
    ASSERT_BETWEEN( 0, local_e, m_ - 1 );
    return local_e + edge_offset_;
  }

  // Returns the id of the local edge (u, v).  The edge must exist.
  EdgeId edge_id(const VertexId &u, const VertexId &v) const
  {
    ASSERT_TRUE( is_vertex_known(u) );
    VertexId local_u = u - vertex_offset_;
    // Perform a binary search for v on the edge heads outgoing from u.
    VertexId *pos = std::lower_bound(
        edge_heads_ + (bucket_begin_offsets_[local_u] - edge_offset_),
        edge_heads_ + (bucket_begin_offsets_[local_u + 1] - edge_offset_),
        v);
    ASSERT_EQ( *pos, v );
    EdgeId e = pos - edge_heads_ + edge_offset_;
    ASSERT_BETWEEN( edge_offset_, e, edge_offset_ + m_ - 1 );
    return e;
  }

  // Returns the pair representation of the edge with the given identifiers.
  // Edge edge(const EdgeId &e) const;  // XXX

  // Returns the global identifier of the edge's tail.
  VertexId edge_tail(const EdgeId &e) const
  { 
    EdgeId local_e = e - edge_offset_;
    ASSERT_BETWEEN( 0, local_e, m_ - 1 );
    return edge_tails_[local_e];
  }

  // Returns the identifier of the edge's head.
  VertexId edge_head(const EdgeId &e) const
  { 
    EdgeId local_e = e - edge_offset_;
    ASSERT_BETWEEN( 0, local_e, m_ - 1 );
    return edge_heads_[local_e];
  }

  // Returns true iff we have information about the star around v in this
  // graph.
  bool is_vertex_known(const VertexId &v) const
  { return (v >= vertex_offset_ and v < vertex_offset_ + n_); }

  // Returns true iff we have information about the edge around v in this
  // graph.
  bool is_edge_known(const EdgeId &e) const
  { return (e >= edge_offset_ and e < edge_offset_ + m_); }

  // Returns the largest local out degree of this graph.
  VertexId max_local_out_degree() const
  { return max_local_out_degree_; }

  // The vertex offset.
  VertexId vertex_offset() const
  { return vertex_offset_; }

  // The edge offset.
  EdgeId edge_offset() const
  { return edge_offset_; }

private:
  // Offset of the vertex ids.
  VertexId vertex_offset_;
  // Offset of the edge ids.
  EdgeId edge_offset_;

  // Number of vertices.
  VertexId n_;
  // Number of edges.
  EdgeId m_;

  // Array with offsets into the edge_{targets, sources}_ array.
  EdgeId *bucket_begin_offsets_;
  // Array that stores the edge tails.
  VertexId *edge_tails_;
  // Array that stores the edge heads.
  VertexId *edge_heads_;

  // Compute and return largest local out degree.
  long ComputeMaxLocalOutDegree() const;

  // The largest local out degree.
  long max_local_out_degree_;

  DISALLOW_COPY_AND_ASSIGN(StaticLocalGraph);

  // For simplicity and performance, we allow friend-access for I/O routines.
  friend int SaveGraphCentrally(
      const GraphFileFormat format,
      const char *filename,
      const StaticLocalGraph *graph,
      const StaticLocalGraphProperties *graph_properties);

  friend class UnsafeStaticLocalGraphAccess;
};


// This class allows the UNSAFE access of the internal structure of
// StaticLocalGraph.
//
// Note that it is your responsibility not to change any of the data or to
// live with undefined behaviour.
class UnsafeStaticLocalGraphAccess
{
public:
  // Initialize the object with a pointer to the given graph.
  UnsafeStaticLocalGraphAccess(const StaticLocalGraph *g) :
    graph_(g)
  {}

  // Return a non-const version of the bucket begin offsets.
  EdgeId *UNSAFE_bucket_begin_offsets()
  { return const_cast<EdgeId*>(graph_->bucket_begin_offsets_); }

  // Return a non-const version of the edge heads.
  VertexId *UNSAFE_edge_heads()
  { return const_cast<VertexId*>(graph_->edge_heads_); }

private:
  // A pointer to the static local graph.
  const StaticLocalGraph *graph_;

  DISALLOW_COPY_AND_ASSIGN(UnsafeStaticLocalGraphAccess);
};


}  // namespace distributed_graph


#endif /* ifndef STATIC_LOCAL_GRAPH_H */

