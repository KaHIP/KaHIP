/* File:   graph_export.h
 * Author: Manuel Holtgrewe <holtgrewe@ira.uka.de>
 */
#ifndef GRAPH_EXPORT_H
#define GRAPH_EXPORT_H

#include "error_codes.h"
#include "graph_types.h"


namespace distributed_graph {


// Forward-Declarations.
//

class StaticLocalGraph;


// Interface
//

// Export the given local graph into the .dot format.
//
// Returns kOk on success and kError on failures.
int ExportGraphCentrallyToDot(const char *filename,
                                  const StaticLocalGraph *graph);


}  // namespace distributed_graph


#endif /* ifndef GRAPH_EXPORT_H */

