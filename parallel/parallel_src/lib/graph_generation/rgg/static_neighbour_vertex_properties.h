/******************************************************************************
 * static_neighbour_vertex_properties.h
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

#ifndef STATIC_NEIGHBOUR_VERTEX_PROPERTIES_H
#define STATIC_NEIGHBOUR_VERTEX_PROPERTIES_H

#include <tr1/unordered_map>

#include <mpi.h>

#include "macros_assertions.h"
#include "macros_common.h"
#include "graph_types.h"

namespace distributed_graph {

// The properties we need to know about a remote vertex.
struct NeighbourVertexInfo
{
  // The number of outgoing edges.
  EdgeCount out_degree;

  // The weight of the vertex.
  VertexWeight weight;

  // The number of internal edges.
  EdgeCount internal_edges;

  // The sum of interface edge weights.
  RggEdgeWeight interface_edge_weight;
};

// A sparse mapping from global vertex identifiers to the vertex weight,
// internal edges and out degree of vertices.
class StaticNeighbourVertexProperties
{
public:
  // Default constructor.
  StaticNeighbourVertexProperties()
  {}

  // Initialize the vertex properties from the given receive buffer.
  //
  // The buffer length is the number of integers, the number of 5-tuples is
  // num / 5.
  //
  // The buffer must contain 5-tuples (v, degree, weight, internal edges,
  // interface edge weight).
  void InitFromReceiveBuffer(const long *buffer, const long num)
  {
    //printf("InitFromReceiveBuffer\n");
    for (int i = 0; i < num; i += 5) {
      NeighbourVertexInfo info;
      info.out_degree = buffer[i + 1];
      info.weight = buffer[i + 2];
      info.internal_edges = buffer[i + 3];
      info.interface_edge_weight = buffer[i + 4];
      infos_[buffer[i]] = info;
    }
  }

  // Returns the NeighbourVertexInfo object.
  const NeighbourVertexInfo &neighbour_info(const VertexId v) const
  {
    ASSERT_TRUE( contains(infos_, v) );
    return infos_.find(v)->second;
  }

  // Return number of outgoing edges for the given global vertex id.
  EdgeCount out_degree(const VertexId v) const
  { return neighbour_info(v).out_degree; }

  // Return weight for the given global vertex id.
  VertexWeight vertex_weight(const VertexId v) const
  { return neighbour_info(v).weight; }

  // Return number of internal edges for the given global vertex id.
  EdgeCount internal_edges(const VertexId v) const
  { return neighbour_info(v).internal_edges; }

  EdgeCount interface_edge_weight(const VertexId v) const
  { return neighbour_info(v).interface_edge_weight; }

private:
  // The info about the neighbour vertex.
  std::tr1::unordered_map<VertexId, NeighbourVertexInfo> infos_;

  DISALLOW_COPY_AND_ASSIGN(StaticNeighbourVertexProperties);
};

}  // namespace distributed_graph

#endif /* ifndef STATIC_NEIGHBOUR_VERTEX_PROPERTIES_H */

