/******************************************************************************
 * graph_export.cpp
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

#include "graph_export.h"  // This file's header.

#include <cstdio> // fprintf()

#include "error_codes.h"
#include "static_local_graph.h"


namespace distributed_graph {


// Implementation.
//


int ExportGraphCentrallyToDot(const char *filename,
                                  const StaticLocalGraph *graph)
{
  // Open file.
  FILE *out = fopen(filename, "w");

  fprintf(out, "graph {\n");
  for (VertexId u = 0, uend = graph->n(); u != uend; ++u) {
    fprintf(out, "%lu; /* degree: %lu */\n", u, graph->out_degree(u));
    for (VertexId i = 0, iend = graph->out_degree(u); i != iend; ++i) {
      VertexId v = graph->neighbour(u, i);
      if (u > v) continue;  // Export each undirected edge once.
      fprintf(out, "%lu -- %lu;\n", u, v);
    }
  }
  fprintf(out, "}\n");

  // Close file again.
  fclose(out);

  return kOk;
}


}

