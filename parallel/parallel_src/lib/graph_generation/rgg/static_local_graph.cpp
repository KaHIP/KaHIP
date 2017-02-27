/******************************************************************************
 * static_local_graph.cpp
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

#include "static_local_graph.h"  // This file's header.

#include <algorithm>  // std::max()
#include <cstdlib>
#include <map>

#include <stdint.h>  // NULL

#include "macros_assertions.h"
#include "util_memory.h"


namespace distributed_graph {


StaticLocalGraph::StaticLocalGraph() :
  vertex_offset_(0), edge_offset_(0), n_(0), m_(0),
  bucket_begin_offsets_(NULL), edge_tails_(NULL), edge_heads_(NULL),
  max_local_out_degree_(-1)
{}


StaticLocalGraph::~StaticLocalGraph()
{
  DeleteArrayIfAssigned(bucket_begin_offsets_);
  DeleteArrayIfAssigned(edge_tails_);
  DeleteArrayIfAssigned(edge_heads_);
}


void StaticLocalGraph::Init(const VertexId vertex_offset,
                            const EdgeId edge_offset,
                            const VertexId n,
                            const EdgeId m,
                            const Edge *edges)
{
  vertex_offset_ = vertex_offset;
  edge_offset_ = edge_offset;
  n_ = n;
  m_ = m;
  max_local_out_degree_ = -1;
  // Allocate memory.
  bucket_begin_offsets_ = new VertexId[n_ + 1];
  edge_tails_ = new VertexId[m_];
  edge_heads_ = new VertexId[m_];

  // Generate adjacency array structure from edge list.
  long i = 0;
  long j = 0;
  bucket_begin_offsets_[0] = edge_offset_;
  while (i < m) {
    // Handle vertices without outgoing edges.
    while (edges[i].tail - vertex_offset > j) {
      ++j;
      bucket_begin_offsets_[j] = i + edge_offset_;

      EdgeCount out_degree = bucket_begin_offsets_[j] -
        bucket_begin_offsets_[j - 1];
      if (out_degree > max_local_out_degree_)
        max_local_out_degree_ = out_degree;
    }

    ASSERT_EQ( edges[i].tail - vertex_offset, j );

    edge_tails_[i] = edges[i].tail;
    edge_heads_[i] = edges[i].head;
    ++i;
  }
  // Handle vertices at the end without outgoing edges.
  while (j < n_) {
    bucket_begin_offsets_[++j] = i + edge_offset_;

    EdgeCount out_degree = bucket_begin_offsets_[j] -
      bucket_begin_offsets_[j - 1];
    if (out_degree > max_local_out_degree_)
      max_local_out_degree_ = out_degree;
  }
  ASSERT_EQ( j, n_ );
  ASSERT_EQ( bucket_begin_offsets_[n_], edge_offset_ + m_ );

  // Finally, compute the maximal out degree.
  max_local_out_degree_ = ComputeMaxLocalOutDegree();
}


void StaticLocalGraph::Init(const VertexId vertex_offset,
                            const EdgeId edge_offset,
                            const VertexId n,
                            const EdgeId m,
                            EdgeId *bucket_begin_offsets,
                            VertexId *edge_tails,
                            VertexId *edge_heads)
{
  vertex_offset_ = vertex_offset;
  edge_offset_ = edge_offset;
  n_ = n;
  m_ = m;
  bucket_begin_offsets_ = bucket_begin_offsets;
  edge_tails_ = edge_tails;
  edge_heads_ = edge_heads;

  /* FIXME(manuel): This does not work if we use it for the gap graph.
  ASSERT_EQ( bucket_begin_offsets_[0], edge_offset );
  ASSERT_EQ( bucket_begin_offsets_[n_], edge_offset + m_ );
  */

  max_local_out_degree_ = ComputeMaxLocalOutDegree();
}


void StaticLocalGraph::Dump(FILE *out) const
{
  fprintf(out, ",-- StaticLocalGraph\n");
  for (VertexId ul = 0; ul < n_; ++ul) {
    VertexId u = ul + vertex_offset_;
    fprintf(out, "| %lu: ", u);
    for (EdgeCount i = 0; i < out_degree(u); ++i) {
      fprintf(out, " %lu", neighbour(u, i));
    }
    fprintf(out, "\n");
  }
  fprintf(out, "`---\n");
}


long StaticLocalGraph::ComputeMaxLocalOutDegree() const
{
  long result = -1;

  for (long i = 1; i < n_ + 1; ++i)
    result = std::max(result,
                      bucket_begin_offsets_[i] - bucket_begin_offsets_[i - 1]);

  return result;
}


void StaticLocalGraph::DumpDegreeHistogram(FILE *out) const
{
  std::map<EdgeCount, long> degree2count;
  for (VertexId ul = 0, ulend = n_; ul < ulend; ++ul) {
    VertexId u = ul + vertex_offset();
    degree2count[out_degree(u)] += 1;
  }
  for (std::map<EdgeCount, long>::const_iterator
         it = degree2count.begin(), itend = degree2count.end();
         it != itend; ++it) {
    fprintf(out, "%lu: %lu  ", it->first, it->second);
  }
}


}  // namespace distributed_graph

