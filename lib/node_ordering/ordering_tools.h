/*
 * Author: Wolfgang Ost
 */

#ifndef ORDERING_TOOLS_H
#define ORDERING_TOOLS_H

#include <iosfwd>
#include <vector>

#include "definitions.h"

// Print the ordering given by 'labels' to the stream 'out' in the format used by scotch:
// the first line contains the number of nodes, each of the following lines contains
// the node index and the ordering.
void print_ordering(std::ostream &out, const std::vector<NodeID> &labels);

class graph_access;
// Compute the adjusted degree for the given node
NodeWeight compute_reachable_set_size(graph_access &graph, NodeID node);

// Compute the number of fill-edges of the given ordering
Count compute_fill(graph_access &graph, const std::vector<NodeID> &ordering);

#endif /* ORDERING_TOOLS_H */
