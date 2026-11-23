/******************************************************************************
 * parallel_graph_io_huge.cpp
 * *
 * Transparent huge pages for the first touch of the graph reader's large arrays.
 * Listed last among the library sources, so the code before it keeps its place.
 *****************************************************************************/

#include <sys/mman.h>
#include <stdint.h>

#include "parallel_graph_io_huge.h"
#include "data_structure/parallel_graph_access.h"

void io_huge_pages::advise(void * p, size_t bytes) {
#ifdef MADV_HUGEPAGE
        const uintptr_t huge  = (uintptr_t)1 << 21;
        const uintptr_t begin = ((uintptr_t)p + huge - 1) & ~(huge - 1);
        const uintptr_t end   = ((uintptr_t)p + bytes) & ~(huge - 1);
        if( end > begin ) madvise((void*)begin, end - begin, MADV_HUGEPAGE); // a hint; THP = never makes it a no-op
#else
        (void)p; (void)bytes;
#endif
}

void io_huge_pages::advise_edges(parallel_graph_access & G, EdgeID entries) {
        if( entries == 0 || entries > G.m_edges.size() ) return;
        advise((void*)G.m_edges.data(), (size_t)entries * sizeof(Edge));
}

ULONG io_huge_pages::entries_lower_bound(const std::vector< ULONG > & all_counts, size_t chunks, ULONG from, ULONG to,
                                         const char * first_chunk_begin, const char * first_line) {
        ULONG lines = 0, entries = 0;
        for( size_t c = 0; c < chunks; c++) {
                const ULONG cl = all_counts[2*c];
                if( cl > 0 && lines >= from && lines + cl <= to + 1 ) entries += all_counts[2*c+1];
                lines += cl;
        }
        // numbers that start before the first line start of the PE's first chunk (same rule as the count scan:
        // a digit after a non-digit; the byte before the chunk is a non-digit or the number started earlier)
        ULONG before = 0;
        for( const unsigned char * q = (const unsigned char *) first_chunk_begin; q < (const unsigned char *) first_line; q++) {
                const bool digit = (unsigned char)(*q - '0') < 10;
                const bool prev_digit = q > (const unsigned char *) first_chunk_begin && (unsigned char)(q[-1] - '0') < 10;
                if( digit && !prev_digit ) before++;
        }
        return entries > before ? entries - before : 0;
}
