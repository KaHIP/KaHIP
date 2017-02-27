/******************************************************************************
 * static_local_graph_properties.cpp
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

#include "static_local_graph_properties.h"  // This file's header.

#include <cstdlib>  // NULL


namespace distributed_graph {


StaticLocalGraphProperties::StaticLocalGraphProperties() :
  n_(0), m_(0), vertex_weights_(NULL), internal_edges_(NULL),
  edge_weights_(NULL), coordinates_(NULL)
{}


StaticLocalGraphProperties::~StaticLocalGraphProperties()
{
  if (vertex_weights_ != NULL)
    delete[] vertex_weights_;
  if (internal_edges_ != NULL)
    delete[] internal_edges_;
  if (edge_weights_ != NULL)
    delete[] edge_weights_;
  if (coordinates_ != NULL)
    delete[] coordinates_;
}


void StaticLocalGraphProperties::InitEmpty(const VertexId vertex_offset,
                                           const EdgeId edge_offset,
                                           const long n,
                                           const long m)
{
  vertex_offset_ = vertex_offset;
  edge_offset_ = edge_offset;
  n_ = n;
  m_ = m;
  vertex_weights_ = new VertexId[n_];
  internal_edges_ = new EdgeCount[n_];
  for (VertexId u = 0; u != n_; ++u) {
    vertex_weights_[u] = 1;
    internal_edges_[u] = 0;
  }
  edge_weights_ = new RggEdgeWeight[m_];
  for (EdgeId e = 0; e != m_; ++e) {
    edge_weights_[e] = 1;
  }
  coordinates_ = NULL;
}


void StaticLocalGraphProperties::InitWithValues(
    const VertexId vertex_offset,
    const EdgeId edge_offset,
    const long n,
    const long m,
    VertexWeight *vertex_weights,
    EdgeCount *internal_edges,
    RggEdgeWeight *edge_weights)
{
  vertex_offset_ = vertex_offset;
  edge_offset_ = edge_offset;
  n_ = n;
  m_ = m;
  vertex_weights_ = vertex_weights;
  internal_edges_ = internal_edges;
  edge_weights_ = edge_weights;
}


void StaticLocalGraphProperties::InitCoordinates(VertexCoordinate *coordinates)
{
  if (coordinates_ != NULL)
    delete[] coordinates_;
  coordinates_ = coordinates;
}


}  // namespace distributed_graph

