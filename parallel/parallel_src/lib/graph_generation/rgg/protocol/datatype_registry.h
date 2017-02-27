// File:   datatype_registry.h
// Author: Manuel Holtgrewe <holtgrewe@ira.uka.de>
//
#ifndef DATATYPE_REGISTRY_H
#define DATATYPE_REGISTRY_H

#include <mpi.h>

namespace distributed_graph {

// This struct is a registry for custom ("derived") MPI datatypes in its static
// members.
//
// Use Register()/Unregister() to register/unregister the datatypes with MPI.
class DatatypeRegistry
{
public:
  // Register all custom datatypes.
  static void Register();

  // Unregister all custom datatypes.
  static void Unregister();

  static MPI::Datatype VERTEX_WITH_COORDINATES;
  static MPI::Datatype IDENTIFIER_PAIR;
  static MPI::Datatype INTEGER_TRIPLE;
  
  static MPI::Op INTEGER_TRIPLE_MIN;

private:
  static void RegisterIdentifierPair();
  static void RegisterIntegerTriple();
};

}  // namespace distributed_graph

#endif // ifndef DATATYPE_REGISTRY_H

