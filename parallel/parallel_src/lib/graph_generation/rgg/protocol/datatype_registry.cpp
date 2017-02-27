// File:   registry.cpp
// Author: Manuel Holtgrewe <holtgrewe@ira.uka.de>
//

#include "datatype_registry.h"  // This file's header.

#include <utility>  // pair<>
#include <mpi.h>

#include "vertex_with_coordinates.h"

namespace distributed_graph {

// Private Interface
//

void MinForTriples(const void *vinvec, void *vinoutvec, int len,
                   const MPI::Datatype& datatype);

// Implementation
//

// Allocate memory for the custom MPI data types and operations.
/*static*/ MPI::Datatype DatatypeRegistry::VERTEX_WITH_COORDINATES = 0;
/*static*/ MPI::Datatype DatatypeRegistry::IDENTIFIER_PAIR = 0;
/*static*/ MPI::Datatype DatatypeRegistry::INTEGER_TRIPLE = 0;
/*static*/ MPI::Op DatatypeRegistry::INTEGER_TRIPLE_MIN = 0;

/*static*/ void DatatypeRegistry::Register()
{
  VertexWithCoordinatesRegisterWithMpi(&VERTEX_WITH_COORDINATES);
  RegisterIdentifierPair();
  RegisterIntegerTriple();
}

/*static*/ void DatatypeRegistry::Unregister()
{
  INTEGER_TRIPLE.Free();
  INTEGER_TRIPLE_MIN.Free();
  IDENTIFIER_PAIR.Free();
  VERTEX_WITH_COORDINATES.Free();
}

/*static*/ void DatatypeRegistry::RegisterIdentifierPair()
{
  IDENTIFIER_PAIR = MPI::INTEGER.Create_contiguous(2);
  IDENTIFIER_PAIR.Commit();
}

/*static*/ void DatatypeRegistry::RegisterIntegerTriple()
{
  INTEGER_TRIPLE = MPI::INTEGER.Create_contiguous(3);
  INTEGER_TRIPLE.Commit();
  INTEGER_TRIPLE_MIN.Init(MinForTriples, true);
}

void MinForTriples(const void *vinvec, void *vinoutvec, int len,
                   const MPI::Datatype& datatype)
{
  // TODO: Only allow INTEGER_TRIPLE!
  const int *invec = static_cast<const int*>(vinvec);
  int *inoutvec = static_cast<int*>(vinoutvec);

  for (int i = 0; i < len; ++i) {
    if (invec[0] < inoutvec[0]) {
      std::copy(invec, invec + 3, inoutvec);
    } else if (invec[0] == inoutvec[0] and invec[1] < inoutvec[1]) {
      std::copy(invec, invec + 3, inoutvec);
    } else if (invec[0] == inoutvec[0] and invec[1] == inoutvec[1] and
               invec[2] < inoutvec[2]) {
      std::copy(invec, invec + 3, inoutvec);
    } else {
      // no copy necessary
    }

    // next
    invec += 3;
    inoutvec += 3;
  }
}

}  // namespace distributed_graph

