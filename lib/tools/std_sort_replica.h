/******************************************************************************
 * std_sort_replica.h
 *
 * std_sort_replica(first, last, comp) produces exactly the same arrangement as
 * libstdc++'s std::sort(first, last, comp) (introsort with median-of-three
 * pivots, unguarded Hoare partition, heap sort fallback and final insertion
 * sort), including the order of elements that compare equal.
 *
 * The only difference is how the Hoare partition finds the elements to swap:
 * the scalar loop pairs the k-th element from the left that does not belong
 * left of the pivot with the k-th element from the right that does not belong
 * right of it, and swaps them while they have not crossed. Here the candidates
 * are collected blockwise without branches (as in BlockQuicksort) and the same
 * pairs are swapped in the same order. Positions scanned ahead but not swapped
 * keep their values, so the scalar loop finishes from the state it would have
 * reached itself.
 *****************************************************************************/

#ifndef STD_SORT_REPLICA_H
#define STD_SORT_REPLICA_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include "tools/cpu_hints.h"

namespace std_sort_replica_detail {

const long SORT_REPLICA_BLOCK = 128;

template<typename T, typename Comp>
inline void move_median_to_first(T* result, T* a, T* b, T* c, Comp comp) {
        if (comp(*a, *b)) {
                if (comp(*b, *c))      std::iter_swap(result, b);
                else if (comp(*a, *c)) std::iter_swap(result, c);
                else                   std::iter_swap(result, a);
        } else if (comp(*a, *c))       std::iter_swap(result, a);
        else if (comp(*b, *c))         std::iter_swap(result, c);
        else                           std::iter_swap(result, b);
}

// libstdc++'s __unguarded_partition(first, last, pivot) with the candidate search done in blocks
template<typename T, typename Comp>
inline T* unguarded_partition(T* first, T* last, const T* pivot, Comp comp) {
        const long B = SORT_REPLICA_BLOCK;
        if (last - first >= 4 * B) {
                unsigned char offl[SORT_REPLICA_BLOCK], offr[SORT_REPLICA_BLOCK];
                T* ls = first;        // left positions < ls are scanned
                T* rs = last;         // right positions >= rs are scanned
                T* lbase = first; T* rbase = last;
                long nl = 0, sl = 0, nr = 0, sr = 0;
                T* hand_first = first; T* hand_last = last;
                const T & p = *pivot;
                while (true) {
                        if (nl == sl && rs - ls >= B) {
                                lbase = ls; nl = 0; sl = 0;
                                for (long i = 0; i < B; i++) {
                                        offl[nl] = (unsigned char) i;
                                        nl += !comp(lbase[i], p);
                                }
                                ls += B;
                        }
                        if (nr == sr && rs - ls >= B) {
                                rbase = rs; nr = 0; sr = 0;
                                for (long i = 1; i <= B; i++) {
                                        offr[nr] = (unsigned char) (i - 1);
                                        nr += !comp(p, rbase[-i]);
                                }
                                rs -= B;
                        }
                        if (nl == sl || nr == sr) break;
                        long n = std::min(nl - sl, nr - sr);
                        for (long q = 0; q < n; q++) {
                                std::iter_swap(lbase + offl[sl + q], rbase - 1 - offr[sr + q]);
                        }
                        hand_first = lbase + offl[sl + n - 1] + 1;
                        hand_last  = rbase - 1 - offr[sr + n - 1];
                        sl += n; sr += n;
                }
                first = hand_first;
                last  = hand_last;
        }
        while (true) {
                while (comp(*first, *pivot)) ++first;
                --last;
                while (comp(*pivot, *last)) --last;
                if (!(first < last)) return first;
                std::iter_swap(first, last);
                ++first;
        }
}

template<typename T, typename Comp>
void introsort_loop(T* first, T* last, std::ptrdiff_t depth_limit, Comp comp) {
        while (last - first > 16) {
                if (depth_limit == 0) {
                        std::partial_sort(first, last, last, comp);
                        return;
                }
                --depth_limit;
                T* mid = first + (last - first) / 2;
                move_median_to_first(first, first + 1, mid, last - 1, comp);
                T* cut = unguarded_partition(first + 1, last, first, comp);
                introsort_loop(cut, last, depth_limit, comp);
                last = cut;
        }
}

template<typename T, typename Comp>
inline void unguarded_linear_insert(T* last, Comp comp) {
        T val = std::move(*last);
        T* next = last - 1;
        while (comp(val, *next)) {
                *last = std::move(*next);
                last = next;
                --next;
        }
        *last = std::move(val);
}

template<typename T, typename Comp>
void insertion_sort(T* first, T* last, Comp comp) {
        if (first == last) return;
        for (T* i = first + 1; i != last; ++i) {
                if (comp(*i, *first)) {
                        T val = std::move(*i);
                        std::move_backward(first, i, i + 1);
                        *first = std::move(val);
                } else {
                        unguarded_linear_insert(i, comp);
                }
        }
}

}

template<typename T, typename Comp>
void std_sort_replica(T* first, T* last, Comp comp) {
        using namespace std_sort_replica_detail;
        if (first == last) return;
        std::ptrdiff_t n = last - first;
        std::ptrdiff_t lg = kahip_floor_log2((uint64_t) n);
        introsort_loop(first, last, lg * 2, comp);
        if (n > 16) {
                insertion_sort(first, first + 16, comp);
                for (T* i = first + 16; i != last; ++i) unguarded_linear_insert(i, comp);
        } else {
                insertion_sort(first, last, comp);
        }
}

#endif
