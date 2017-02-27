/******************************************************************************
 * graph_io.cpp
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

#include "graph_io.h"  // This file's header.

#include <cmath>    // ceil()
#include <cstdio>
#include <cstring>
#include <memory>
#include <vector>
#include <fstream>
#include <iostream>

#include <mpi.h>

#include "error_codes.h"
//#include "instrumentation.h"
#include "macros_assertions.h"
#include "macros_common.h"
#include "static_distributed_graph.h"
#include "static_local_graph.h"
#include "static_local_graph_properties.h"
#include "util.h"

// The maximal line length supported in graph files.
#define BUFFER_LEN (128*1024)

namespace distributed_graph {

// Local Definitions.
//

// Magic number of the header.
const int kDaaMagicNumber[3] = {'D', 'A', 'A'};
// Number of int values in the header.
const int kDaaHeaderIntCount = 7;
// The current file version.
const int kDaaCurrentVersion = 2;

// METIS Graph format value:  No weights.
const int kDaaFormatNoWeights = 0;
// METIS Graph format value:  Edge weights only.
const int kDaaFormatEdgeWeights = 1;
// METIS Graph format value:  Vertex weights only.
const int kDaaFormatVertexWeights = 10;
// METIS Graph format value:  Both vertex and edge weights.
const int kDaaFormatVertexEdgeWeights = 11;

// Magic number of the coordinate header.
const int kCoordsMagicNumber[3] = {'C', 'D', 'S'};
// Number of int values in the coords header.
const int kCoordsHeaderIntCount = 5;
// The current coords file version
const int kCoordsCurrentVersion = 1;


int SaveGraphDistributedlyMetis(const char *filename,
                               const StaticDistributedGraph *graph)
{

  //VertexId n = graph->global_n();
  //EdgeId m = graph->global_m();

  int rank = graph->communicator()->Get_rank();
  int size = graph->communicator()->Get_size();

  const VertexId *vertex_splitters = graph->vertex_splitters();

  // Get a shortcut to the current static local graph.
  const StaticLocalGraph *local_graph = graph->local_graph();
  UnsafeStaticLocalGraphAccess local_access(local_graph);

  if( rank == 0) {
        std::ofstream f(filename);
        f << graph->global_n() <<  " " <<  graph->global_m()/2 << std::endl;
        f.close();
  }

  for( int peID = 0; peID < size; peID++) {
          graph->communicator()->Barrier();
          if( rank == peID ) {
                  std::ofstream f;   
                  f.open(filename, std::ofstream::out | std::ofstream::app);
                  for( long i = 0, iend = local_graph->n(); i < iend; i++) {
                          VertexId u = i+vertex_splitters[peID];
                          for( long j = 0, jend = local_graph->out_degree(u); j < jend; j++) {
                                  //EdgeId e = local_graph->out_edge(u,j);
                                  VertexId v = local_graph->neighbour(u,j);
                                  f     << (v+1) << " ";
                          }
                          f << std::endl;
                  }
                  f.close();
          }
  }

  graph->communicator()->Barrier();

  return kOk;
}


int StoreCoordinatesDistributedlyMetis(const char *filename,
                                     const StaticDistributedGraph *graph,
                                     const VertexCoordinate *coordinates)
{
        const StaticLocalGraph *local_graph = graph->local_graph();
        int rank = graph->communicator()->Get_rank();
        int size = graph->communicator()->Get_size();

        if( rank == 0) {
                std::ofstream f(filename);
                for( long i = 0, iend = local_graph->n(); i < iend; i++) {
                        f << coordinates[i*3 + 0] << " " << coordinates[i*3 + 1] << " " << coordinates[i*3 + 2] << std::endl;
                }
                f.close();
        }


        for( int peID = 1; peID < size; peID++) {
                graph->communicator()->Barrier();
                if( rank == peID ) {
                        std::ofstream f;   
                        f.open(filename, std::ofstream::out | std::ofstream::app);
                        for( long i = 0, iend = local_graph->n(); i < iend; i++) {
                                f << coordinates[i*3 + 0] << " " << coordinates[i*3 + 1] << " " << coordinates[i*3 + 2] << std::endl;
                        }
                        f.close();
                }
        }

        graph->communicator()->Barrier();

        return kOk;
}



}  // namespace distributed_graph

