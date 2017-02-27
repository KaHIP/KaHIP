/* File:   static_local_graph_properties.cpp
 * Author: Manuel Holtgrewe <holtgrewe@ira.uka.de>
 *
 * Implementation of the StaticLocalGraphProperties class.
 */

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

