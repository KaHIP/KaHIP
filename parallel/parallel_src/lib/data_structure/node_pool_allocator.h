/******************************************************************************
 * node_pool_allocator.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef NODE_POOL_ALLOCATOR_PARHIP_H
#define NODE_POOL_ALLOCATOR_PARHIP_H

#include <cstddef>
#include <new>
#include <vector>

// Free list of fixed-size blocks, carved from 64 KiB chunks that are kept until the process ends.
// ParHIP runs one thread per process, so no locking.
template < std::size_t Size >
class node_pool {
public:
        static void * get() {
                node_pool & p = instance();
                if( p.m_free != NULL ) {
                        free_block * b = p.m_free;
                        p.m_free = b->next;
                        return b;
                }
                if( p.m_cur == p.m_end ) p.new_chunk();
                void * r = p.m_cur;
                p.m_cur += block_size();
                return r;
        }

        static void put( void * ptr ) {
                node_pool & p = instance();
                free_block * b = static_cast< free_block * >(ptr);
                b->next  = p.m_free;
                p.m_free = b;
        }

private:
        struct free_block { free_block * next; };

        static std::size_t block_size() {
                std::size_t s = Size < sizeof(free_block) ? sizeof(free_block) : Size;
                const std::size_t a = alignof(std::max_align_t);
                return (s + a - 1) / a * a;
        }

        node_pool() : m_free(NULL), m_cur(NULL), m_end(NULL) {}

        static node_pool & instance() {
                static node_pool p;
                return p;
        }

        void new_chunk() {
                const std::size_t count = 65536 / block_size() + 1;
                char * c = static_cast< char * >(::operator new(count * block_size()));
                m_cur = c;
                m_end = c + count * block_size();
        }

        free_block * m_free;
        char * m_cur;
        char * m_end;
};

// Stateless allocator: single objects (the nodes of node-based containers) come from node_pool,
// arrays (bucket arrays) from operator new. The allocator does not influence the iteration order
// of std::unordered_map, which depends only on the hash values and the insertion sequence.
template < typename T >
struct node_pool_allocator {
        typedef T value_type;

        node_pool_allocator() {}
        template < typename U > node_pool_allocator( const node_pool_allocator< U > & ) {}

        T * allocate( std::size_t n ) {
                if( n == 1 ) return static_cast< T * >(node_pool< sizeof(T) >::get());
                return static_cast< T * >(::operator new(n * sizeof(T)));
        }

        void deallocate( T * p, std::size_t n ) {
                if( n == 1 ) node_pool< sizeof(T) >::put(p);
                else         ::operator delete(p);
        }
};

template < typename T, typename U >
bool operator==( const node_pool_allocator< T > &, const node_pool_allocator< U > & ) { return true; }
template < typename T, typename U >
bool operator!=( const node_pool_allocator< T > &, const node_pool_allocator< U > & ) { return false; }

#endif
