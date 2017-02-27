/******************************************************************************
 * graph_export.h
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

