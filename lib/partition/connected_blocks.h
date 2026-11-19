/******************************************************************************
 * connected_blocks.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef CONNECTED_BLOCKS_4R8T2W6Q
#define CONNECTED_BLOCKS_4R8T2W6Q

#include "data_structure/graph_access.h"
#include "partition/partition_config.h"

// Returns true if G is connected (an empty graph counts as connected).
bool graph_is_connected(graph_access & G);

// Post-processing for connected blocks (used by kaffpa --connected_blocks and
// the kaffpa_connected_blocks library call). Expects G to be partitioned with
// partition_config.connected_blocks set. Moves every component of a block except
// the heaviest one to the lightest adjacent block, greedily rebalances by moving
// boundary nodes that do not disconnect their block, and finishes with a k-way
// refinement pass that keeps blocks connected. Each block is connected afterwards
// if the input graph is connected.
void enforce_connected_blocks(PartitionConfig & partition_config, graph_access & G);

#endif /* end of include guard: CONNECTED_BLOCKS_4R8T2W6Q */
