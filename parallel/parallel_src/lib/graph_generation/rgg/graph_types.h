/******************************************************************************
 * graph_types.h
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


#ifndef GRAPH_TYPES_H
#define GRAPH_TYPES_H

#include <stdint.h>  // int

#include <mpi.h>


namespace distributed_graph {


// Enum specifying a graph file format.
enum GraphFileFormat
{
  GRAPH_FORMAT_INVALID,
  GRAPH_FORMAT_METIS,   // METIS text format.
  GRAPH_FORMAT_DAA,     // Binary adjacency array.
  GRAPH_FORMAT_DDSG     // DDSG graph format.
};


const int kInvalidProcessId = -1;
// Constant with the value of the invalid vertex identifier.
const int kInvalidVertexId = -1;
// Constant with the value of the invalid edge identifier.
const int kInvalidEdgeId = -1;
// Constant for the key parameter in MPI::Intracomm::Split() for "ignore".
const int kMpiKeyIgnore = -1;


// Define a generic type we use for all kinds of identifiers below.  We will
// also use this type for weights and counts.
typedef long IdentifierType;

typedef IdentifierType VertexId;
typedef IdentifierType EdgeId;
typedef IdentifierType PartitionId;
typedef IdentifierType ProcessId;

// Define a type for counting things.
typedef IdentifierType CountType;

typedef CountType VertexWeight;
typedef CountType VertexCount;
typedef CountType RggEdgeWeight;
typedef CountType EdgeCount;

// We use a floating point value for rating edges by rating.
// TODO(manuel): If we get overflows then use a double.
// TODO(manuel): Rename to EdgeRating.
typedef float IndexEdgeRating;


// An entry in the coordinate vector of a vertice's coordinate.
typedef float VertexCoordinate;


// Edge representation as a pair of vertex identifiers.
//
// This is more or less a copy of std::pair<,>.  However, we feel that having
// members called tail/head instead of first/second helps clearity of the code.
struct Edge
{
  // Initialize tail and head to 0.
  Edge() :
    tail(kInvalidVertexId), head(kInvalidVertexId)
  {}

  // Initialize tail and head with the given values.
  Edge(const VertexId &t, const VertexId &h) :
    tail(t), head(h)
  {}

  // Test two edges for equality.
  bool operator==(const Edge &right) const
  { return tail == right.tail and head == right.head; }

  // Test two edges for inequality.
  bool operator!=(const Edge &right) const
  { return tail != right.tail or head != right.head; }

  // Strict Lesser Than Ordering;
  bool operator<(const Edge &right) const
  { return (tail < right.tail) or (tail == right.tail and head < right.head); }

  VertexId tail;
  VertexId head;
};


// Extract the tail from an edge.
inline
VertexId ExtractEdgeTail(const Edge &edge)
{ return edge.tail; }

// Compare two edges by tail vertex id.
inline
bool CompareEdgesByTail(const Edge &left, const Edge &right)
{ return left.tail < right.tail; }

// Compare two edges lexicographically.
inline 
bool CompareEdgesLexicographically(const Edge &left, const Edge &right)
{
  return (left.tail < right.tail) or
    (left.tail == right.tail and left.head < right.head);
}


// Register the Edge type with MPI.
//
// The type will also be committed.
//
// TODO(manuel): Maybe we should have an extra header for our MPI types?
inline
void CreateMpiStructForEdge(MPI::Datatype *newtype)
{
  // Setup block lengths and types.
  int blens[] = {1, 1};
  MPI::Datatype datatypes[] = {MPI::LONG, MPI::LONG};
  // Compute displacements.
  Edge e;
  MPI::Aint address_e = MPI::Get_address(&e);
  MPI::Aint address_e_tail = MPI::Get_address(&(e.tail));
  MPI::Aint address_e_head = MPI::Get_address(&(e.head));
  MPI::Aint displacements[2];
  displacements[0] = address_e_tail - address_e;
  displacements[1] = address_e_head - address_e;

  // Create and commit the new type.
  *newtype = MPI::Datatype::Create_struct(2, blens, displacements, datatypes);
  newtype->Commit();
}



}  // namespace distributed_graph


#endif /* ifndef GRAPH_TYPES_H */

