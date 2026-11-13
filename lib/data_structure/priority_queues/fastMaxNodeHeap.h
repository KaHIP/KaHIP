/******************************************************************************
 * fastMaxNodeHeap.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef FAST_MAX_NODE_HEAP_39CK1B8I
#define FAST_MAX_NODE_HEAP_39CK1B8I

#include <limits>
#include <vector>
#include <unordered_map>
#ifndef _WIN32
#include <execinfo.h>
#endif

#include "data_structure/priority_queues/priority_queue_interface.h"
#include "data_structure/node_hash_map.h"
#include "data_structure/stamped_node_map.h"

typedef int Key;

// Binary max-heap of nodes. Each heap entry holds the key and the node, and a node map gives the
// heap position of each node. The sift operations make the same comparisons and leave every entry
// at the same position as the former implementation (entries swapped one level at a time, with an
// element array between heap and node map), so ties are broken identically.
class fastMaxNodeHeap final : public priority_queue_interface {
        public:
                fastMaxNodeHeap() {};
                ~fastMaxNodeHeap() override = default;
                fastMaxNodeHeap(fastMaxNodeHeap &&) = default;
                fastMaxNodeHeap & operator=(fastMaxNodeHeap &&) = default;

                NodeID size() override;
                bool empty() override;

                bool contains(NodeID node) override;
                void insert(NodeID id, Gain gain) override;

                NodeID deleteMax() override;
                void deleteNode(NodeID node) override;
                NodeID maxElement() override;
                Gain maxValue() override;

                void decreaseKey(NodeID node, Gain gain) override;
                void increaseKey(NodeID node, Gain gain) override;
                void changeKey(NodeID node, Gain gain) override;
                Gain getKey(NodeID node) override;

        private:
                struct entry {
                        Key    key;
                        NodeID node;
                };
                std::vector< entry >   m_heap;     // the heap
                stamped_node_map<int>  m_position; // position of each node in m_heap

                void siftUp( int pos );
                void siftDown( int pos );
                void place( int pos, const entry & x ) {
                        m_heap[pos] = x;
                        m_position.present(x.node) = pos;
                }
};

inline Gain fastMaxNodeHeap::maxValue() {
        return m_heap[0].key;
};

inline NodeID fastMaxNodeHeap::maxElement() {
        return m_heap[0].node;
};

// Moves the entry at pos down: while it has two children and one of them has a key >= its key,
// it changes places with the larger child (the right one if both are equal); with only a left child
// it changes places if that key is larger.
inline void fastMaxNodeHeap::siftDown( int pos ) {
        const entry x    = m_heap[pos];
        const int   size = m_heap.size();
        while(true) {
                const int lhsChild = 2*pos+1;
                const int rhsChild = 2*pos+2;
                int next;
                if( rhsChild < size ) {
                        const Key lhsKey = m_heap[lhsChild].key;
                        const Key rhsKey = m_heap[rhsChild].key;
                        if( lhsKey < x.key && rhsKey < x.key) break;
                        next = lhsKey > rhsKey ? lhsChild : rhsChild;
                } else if( lhsChild < size && x.key < m_heap[lhsChild].key ) {
                        next = lhsChild;
                } else break;
                place(pos, m_heap[next]);
                pos = next;
        }
        place(pos, x);
}

// Moves the entry at pos up while its parent has a smaller key.
inline void fastMaxNodeHeap::siftUp( int pos ) {
        const entry x = m_heap[pos];
        while( pos > 0 ) {
                const int parentPos = (pos-1)/2;
                if( !(m_heap[parentPos].key < x.key) ) break;
                place(pos, m_heap[parentPos]);
                pos = parentPos;
        }
        place(pos, x);
}

inline NodeID fastMaxNodeHeap::size() {
        return m_heap.size();
}

inline bool fastMaxNodeHeap::empty( ) {
        return m_heap.empty();
}

inline void fastMaxNodeHeap::insert(NodeID node, Gain gain) {
        if( m_position.find(node) == m_position.end() ) {
                int heap_size = m_heap.size();
                m_heap.push_back( entry{gain, node} );
                m_position[node] = heap_size;
                siftUp( heap_size );
        }
}

inline void fastMaxNodeHeap::deleteNode(NodeID node) {
        int heap_index = m_position.present(node);
        m_position.erase(node);

        const int last = m_heap.size() - 1;
        if( heap_index != last ) place(heap_index, m_heap[last]);
        m_heap.pop_back();

        if( m_heap.size() > 1 && heap_index < (int)m_heap.size() ) {
                //fix the max heap property
                siftDown(heap_index);
                siftUp(heap_index);
        }
};

inline NodeID fastMaxNodeHeap::deleteMax() {
        if( m_heap.size() > 0) {
                NodeID node = m_heap[0].node;
                m_position.erase(node);

                const int last = m_heap.size() - 1;
                if( last != 0 ) place(0, m_heap[last]);
                m_heap.pop_back();

                if( m_heap.size() > 1) {
                        //fix the heap property
                        siftDown(0);
                }

                return node;
        }

        return -1;
}

inline void fastMaxNodeHeap::changeKey(NodeID node, Gain gain) {
        Gain old_gain = m_heap[m_position.present(node)].key;
        if( old_gain > gain ) {
                decreaseKey(node, gain);
        } else if ( old_gain < gain ) {
                increaseKey(node, gain);
        }
};

inline void fastMaxNodeHeap::decreaseKey(NodeID node, Gain gain) {
        ASSERT_TRUE(m_position.find(node) != m_position.end());
        int heap_idx = m_position.present(node);
        m_heap[heap_idx].key = gain;
        siftDown(heap_idx);
}

inline void fastMaxNodeHeap::increaseKey(NodeID node, Gain gain) {
        ASSERT_TRUE(m_position.find(node) != m_position.end());
        int heap_idx = m_position.present(node);
        m_heap[heap_idx].key = gain;
        siftUp(heap_idx);
}

inline Gain fastMaxNodeHeap::getKey(NodeID node) {
        return m_heap[m_position.present(node)].key;
};


inline bool fastMaxNodeHeap::contains(NodeID node) {
       return m_position.find(node) != m_position.end();
}

#endif
