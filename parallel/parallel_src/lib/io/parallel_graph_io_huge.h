/******************************************************************************
 * parallel_graph_io_huge.h
 * *
 * Transparent huge pages for the first touch of the graph reader's large arrays
 * (a hint to the kernel only: values, order and results are unchanged).
 *****************************************************************************/

#ifndef PARALLEL_GRAPH_IO_HUGE_H
#define PARALLEL_GRAPH_IO_HUGE_H

#include <cstddef>
#include <vector>

#include "definitions.h"

class parallel_graph_access;

struct io_huge_pages {
        // madvise(MADV_HUGEPAGE) on the 2 MiB-aligned interior of [p, p + bytes); only for ranges that are
        // written in full before they are read, so that no huge page holds bytes that would stay untouched
        static void advise(void * p, size_t bytes);

        // the first entries elements of G's edge array (sized by start_construction, not yet touched)
        static void advise_edges(parallel_graph_access & G, EdgeID entries);

        // a lower bound of the adjacency entries the PE with node range [from, to] reads, from the per-chunk
        // counts of readGraphWeightedFlexible (chunk_lines[c], chunk_entries(c)): the chunks that start a line
        // and whose lines all lie in [from, to], minus the entries before the first line start of the first
        // such chunk (they belong to the previous PE's last line); chunks without a line start are left out
        static ULONG entries_lower_bound(const std::vector< ULONG > & all_counts, size_t chunks, ULONG from, ULONG to,
                                         const char * first_chunk_begin, const char * first_line);
};

#endif
