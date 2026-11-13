/******************************************************************************
 * cpu_hints.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 * Portable wrappers for compiler intrinsics: software prefetch (a hint only,
 * it never changes results) and the floor of log2 of a 64-bit integer.
 * GCC and Clang use their builtins, MSVC its intrinsics; other compilers get
 * a no-op prefetch and a plain loop.
 *****************************************************************************/

#ifndef CPU_HINTS_H
#define CPU_HINTS_H

#include <cstdint>

#if defined(__GNUC__) || defined(__clang__)
#define KAHIP_PREFETCH(addr)       __builtin_prefetch(addr)
#define KAHIP_PREFETCH_WRITE(addr) __builtin_prefetch((addr), 1)
#elif defined(_MSC_VER) && (defined(_M_X64) || defined(_M_IX86))
#include <xmmintrin.h>
#define KAHIP_PREFETCH(addr)       _mm_prefetch((const char*)(addr), _MM_HINT_T0)
#define KAHIP_PREFETCH_WRITE(addr) _mm_prefetch((const char*)(addr), _MM_HINT_T0)
#else
#define KAHIP_PREFETCH(addr)       ((void)0)
#define KAHIP_PREFETCH_WRITE(addr) ((void)0)
#endif

#if defined(_MSC_VER) && !defined(__clang__) && defined(_M_X64)
#include <intrin.h>
#endif

// floor(log2(n)) for n > 0 (the index of the highest set bit)
inline int kahip_floor_log2(uint64_t n) {
#if defined(__GNUC__) || defined(__clang__)
        return 63 - __builtin_clzll((unsigned long long) n);
#elif defined(_MSC_VER) && defined(_M_X64)
        unsigned long idx;
        _BitScanReverse64(&idx, n);
        return (int) idx;
#else
        int r = 0;
        while (n >>= 1) r++;
        return r;
#endif
}

#endif
