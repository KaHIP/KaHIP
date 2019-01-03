/******************************************************************************
 * tabu_moves_queue.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef TABU_MOVES_PQ_EM8YJPA9
#define TABU_MOVES_PQ_EM8YJPA9

#include <limits>
#include <queue>

#include "data_structure/matrix/normal_matrix.h"

struct TabuTimePair {
        int time;
        NodeID node;
        PartitionID block;
};
struct comparePair{
        bool operator() (TabuTimePair& x, TabuTimePair& y)
        { return x.time > y.time; }
};


//maxheap 
typedef std::priority_queue< TabuTimePair, std::vector< TabuTimePair >, comparePair > PQ;

class tabu_moves_queue  {
        public:
                tabu_moves_queue( ); 
                virtual ~tabu_moves_queue() { };

                NodeID size();  
                bool empty();

                void insert(NodeID node, PartitionID block, int time); 
                int minValue();
                std::pair<NodeID, PartitionID> deleteMin();

                bool contains(NodeID node, PartitionID block);
        private:
                PQ m_priority_queue;
};

inline tabu_moves_queue::tabu_moves_queue() {
}

inline NodeID tabu_moves_queue::size() {
        return m_priority_queue.size();  
}

inline void tabu_moves_queue::insert(NodeID node, PartitionID block, int time) {
        TabuTimePair ttp;
        ttp.node  = node;
        ttp.time  = time;
        ttp.block = block;
        m_priority_queue.push(ttp);
}

inline bool tabu_moves_queue::empty( ) {
        return m_priority_queue.empty();        
}

inline Gain tabu_moves_queue::minValue( ) {
        return m_priority_queue.top().time;        
}

inline std::pair<NodeID, PartitionID> tabu_moves_queue::deleteMin() {
        std::pair< NodeID, PartitionID > p;
        p.first  = m_priority_queue.top().node;
        p.second = m_priority_queue.top().block;
        m_priority_queue.pop();
        return p;        
}

#endif /* end of include guard: BUCKET_PQ_EM8YJPA9 */
