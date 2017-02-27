/******************************************************************************
 * datatype_registry.h
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

