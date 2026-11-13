/******************************************************************************
 * small_block_alloc.h
 *
 * Replacement of the global operator new/delete for the single-threaded
 * kaffpa binary: small blocks (up to 512 bytes) come from per-size-class
 * free lists inside one reserved virtual region (class = address range, so no
 * per-block header); larger blocks go to malloc. KaHIP never orders or hashes
 * by pointer value, so the allocator cannot change any result.
 *
 * Include in exactly one translation unit (app/kaffpa.cpp).
 *****************************************************************************/

#ifndef SMALL_BLOCK_ALLOC_H
#define SMALL_BLOCK_ALLOC_H

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <new>
#include <sys/mman.h>

namespace small_block_alloc {

static const size_t GRANULE     = 16;
static const size_t NCLASSES    = 32;            // sizes 16, 32, ..., 512
static const size_t CLASS_SHIFT = 32;            // 4 GiB of address space per class
static const size_t MAX_SMALL   = GRANULE * NCLASSES;

struct free_block { free_block * next; };

struct state_t {
        char *       base;                       // start of the reserved region, nullptr if unavailable
        bool         initialized;
        char *       bump[NCLASSES];
        char *       limit[NCLASSES];
        free_block * free_list[NCLASSES];
};

static state_t g_state;

static inline void init() {
        g_state.initialized = true;
        void * p = mmap(nullptr, NCLASSES << CLASS_SHIFT, PROT_READ | PROT_WRITE,
                        MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
        if(p == MAP_FAILED) {
                g_state.base = nullptr;
                return;
        }
        g_state.base = (char*) p;
        for(size_t c = 0; c < NCLASSES; c++) {
                g_state.bump[c]      = g_state.base + (c << CLASS_SHIFT);
                g_state.limit[c]     = g_state.bump[c] + ((size_t)1 << CLASS_SHIFT);
                g_state.free_list[c] = nullptr;
        }
}

static inline void * allocate(size_t size) {
        if(size <= MAX_SMALL) {
                if(!g_state.initialized) init();
                if(g_state.base) {
                        size_t c = size == 0 ? 0 : (size - 1) / GRANULE;
                        free_block * b = g_state.free_list[c];
                        if(b) {
                                g_state.free_list[c] = b->next;
                                return b;
                        }
                        size_t sz = (c + 1) * GRANULE;
                        if(g_state.bump[c] + sz <= g_state.limit[c]) {
                                void * r = g_state.bump[c];
                                g_state.bump[c] += sz;
                                return r;
                        }
                }
        }
        void * r = std::malloc(size == 0 ? 1 : size);
        if(!r) throw std::bad_alloc();
        return r;
}

static inline void deallocate(void * p) {
        if(!p) return;
        size_t off = (size_t)((char*)p - g_state.base);
        if(g_state.base && (char*)p >= g_state.base && off < (NCLASSES << CLASS_SHIFT)) {
                size_t c = off >> CLASS_SHIFT;
                free_block * b = (free_block*) p;
                b->next = g_state.free_list[c];
                g_state.free_list[c] = b;
                return;
        }
        std::free(p);
}

}

void * operator new(size_t size) { return small_block_alloc::allocate(size); }
void * operator new[](size_t size) { return small_block_alloc::allocate(size); }
void * operator new(size_t size, const std::nothrow_t &) noexcept {
        try { return small_block_alloc::allocate(size); } catch(...) { return nullptr; }
}
void * operator new[](size_t size, const std::nothrow_t &) noexcept {
        try { return small_block_alloc::allocate(size); } catch(...) { return nullptr; }
}
void operator delete(void * p) noexcept { small_block_alloc::deallocate(p); }
void operator delete[](void * p) noexcept { small_block_alloc::deallocate(p); }
void operator delete(void * p, const std::nothrow_t &) noexcept { small_block_alloc::deallocate(p); }
void operator delete[](void * p, const std::nothrow_t &) noexcept { small_block_alloc::deallocate(p); }
void operator delete(void * p, size_t) noexcept { small_block_alloc::deallocate(p); }
void operator delete[](void * p, size_t) noexcept { small_block_alloc::deallocate(p); }

#endif
