/******************************************************************************
 * ordered_edge_map_test.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

// Inserts random edge sequences into a hashed_graph (std::unordered_map) and into an
// ordered_edge_map and checks that both iterate the same edges with the same orientation and
// weight in the same order. The coarse graphs of ParHIP depend on this order.
// Usage: ordered_edge_map_test [number of sequences, default 300]

#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

#include "data_structure/ordered_edge_map.h"

int main( int argc, char ** argv ) {
#ifdef ORDERED_EDGE_MAP_EMULATES_UNORDERED_MAP
        printf("ordered_edge_map: emulation of the unordered_map order\n");
#else
        printf("ordered_edge_map: std::unordered_map\n");
#endif
        int sequences = argc > 1 ? atoi(argv[1]) : 300;
        std::mt19937_64 rng(1);
        int bad = 0;
        for( int t = 0; t < sequences; t++ ) {
                // two thirds small sequences (many rehash epochs relative to the size),
                // one third large ones; k as in the contraction, plus large k
                bool small     = t < sequences * 2 / 3;
                NodeID n       = 1 + rng() % (small ? 200 : 200000);
                size_t inserts = rng() % (small ? 2000 : 600000);
                NodeID k       = (t % 3 == 0) ? n : (rng() % 3 == 0 ? (1ULL << 40) + 7 : n + rng() % 1000);

                hashed_graph hG;
                ordered_edge_map M(k);
                for( size_t i = 0; i < inserts; i++ ) {
                        NodeID s = rng() % n, d = rng() % n;
                        if( s == d ) continue;
                        NodeWeight w = 1 + rng() % 5;
                        hashed_edge he;
                        he.k = k; he.source = s; he.target = d;
                        hG[he].weight += w;
                        M.weight(s, d) += w;
                }

                std::vector< NodeID > a, b;
                for( hashed_graph::iterator it = hG.begin(); it != hG.end(); ++it ) {
                        a.push_back(it->first.source);
                        a.push_back(it->first.target);
                        a.push_back(it->second.weight);
                }
                M.for_each([&](NodeID s, NodeID d, NodeWeight w) {
                        b.push_back(s);
                        b.push_back(d);
                        b.push_back(w);
                });
                if( a != b || M.size() != hG.size() ) {
                        if( bad < 5 ) printf("mismatch in sequence %d (n %llu, inserts %zu)\n",
                                             t, (unsigned long long)n, inserts);
                        bad++;
                }
        }
        printf("%d sequences, %d mismatches\n", sequences, bad);
        return bad == 0 ? 0 : 1;
}
