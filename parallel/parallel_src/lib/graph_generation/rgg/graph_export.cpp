/* File:   graph_export.cpp
 * Author: Manuel Holtgrewe <holtgrewe@ira.uka.de>
 */

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

