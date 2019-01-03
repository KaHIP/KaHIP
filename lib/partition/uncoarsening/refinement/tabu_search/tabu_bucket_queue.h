/******************************************************************************
 * tabu_bucket_queue.h 
 *
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef TABU_BUCKET_PQ_EM8YJPA9
#define TABU_BUCKET_PQ_EM8YJPA9

#include <limits>

//this PQ is specalized for Tabu Search, it only contains non-tabu moves
//there is a second PQ that contains tabu moves
#include "data_structure/matrix/normal_matrix.h"
#include "data_structure/priority_queues/priority_queue_interface.h"
#include "random_functions.h"

class tabu_bucket_queue  {
        public:
                tabu_bucket_queue( PartitionConfig & config, const EdgeWeight & gain_span, NodeID number_of_nodes ); 

                virtual ~tabu_bucket_queue() { delete m_queue_index; delete m_gains;};

                NodeID size();  
                void insert(NodeID id, PartitionID block, Gain gain); 
                bool empty();

                Gain maxValue();
                std::pair<NodeID, PartitionID> maxElement();
                std::pair<NodeID, PartitionID> deleteMax();

                void decreaseKey(NodeID node, PartitionID block, Gain newGain);
                void increaseKey(NodeID node, PartitionID block, Gain newGain);

                void changeKey(NodeID element, PartitionID block, Gain newKey);
                Gain getKey(NodeID element, PartitionID block);
                void deleteNode(NodeID node, PartitionID block);

                bool contains(NodeID node, PartitionID block);
        private:
                normal_matrix* m_queue_index;
                normal_matrix* m_gains;
                NodeID         m_elements;
                EdgeWeight     m_gain_span;
                unsigned       m_max_idx; //points to the non-empty bucket with the largest gain
                
                std::vector< std::vector< std::pair<NodeID, PartitionID> > > m_buckets;
};

inline tabu_bucket_queue::tabu_bucket_queue( PartitionConfig & config, 
                                             const EdgeWeight & gain_span_input, 
                                             NodeID number_of_nodes ) {
        m_elements    = 0;
        m_gain_span   = gain_span_input;
        m_max_idx     = 0;
        m_queue_index = new normal_matrix( number_of_nodes, config.k, NOTINQUEUE);
        m_gains       = new normal_matrix( number_of_nodes, config.k, NOTINQUEUE);
        m_buckets.resize(2*m_gain_span+1);
}

inline NodeID tabu_bucket_queue::size() {
        return m_elements;  
}

inline void tabu_bucket_queue::insert(NodeID node, PartitionID block, Gain gain) {
        unsigned address = gain + m_gain_span;
        if(address > m_max_idx) {
                m_max_idx = address; 
        }
       
        std::pair< NodeID, PartitionID > p;
        p.first  = node;
        p.second = block;

        m_buckets[address].push_back( p ); 
        m_queue_index->set_xy(node, block, m_buckets[address].size() - 1); //store position
        m_gains->set_xy(node, block, gain);
 
        m_elements++;
}

inline bool tabu_bucket_queue::empty( ) {
        return m_elements == 0;        
}

inline Gain tabu_bucket_queue::maxValue( ) {
        return m_max_idx - m_gain_span;        
}

inline std::pair<NodeID, PartitionID> tabu_bucket_queue::maxElement( ) {
        return m_buckets[m_max_idx].back();        
}

inline std::pair<NodeID, PartitionID> tabu_bucket_queue::deleteMax() {
       unsigned rnd_idx = random_functions::nextInt(0, m_buckets[m_max_idx].size()-1);
       swap(m_buckets[m_max_idx][rnd_idx], m_buckets[m_max_idx].back());
       m_queue_index->set_xy(m_buckets[m_max_idx][rnd_idx].first, m_buckets[m_max_idx][rnd_idx].second, rnd_idx); 

       std::pair< NodeID, PartitionID > p;
       p = m_buckets[m_max_idx].back();
       m_buckets[m_max_idx].pop_back();

       m_queue_index->set_xy(p.first, p.second, NOTINQUEUE); //erase(node, block);
       m_gains->set_xy(p.first, p.second, NOTINQUEUE);

       if( m_buckets[m_max_idx].size() == 0 ) {
             //update max_idx
             while( m_max_idx != 0 )  {
                     m_max_idx--;
                     if(m_buckets[m_max_idx].size() > 0) {
                        break;
                     }
             }
       }

       m_elements--;
       return p;        
}

inline void tabu_bucket_queue::decreaseKey(NodeID node, PartitionID block, Gain new_gain) {
        changeKey( node, block, new_gain );
}

inline void tabu_bucket_queue::increaseKey(NodeID node, PartitionID block, Gain new_gain) {
        changeKey( node, block, new_gain );
}

inline void tabu_bucket_queue::changeKey(NodeID node, PartitionID block, Gain new_gain) {
        deleteNode(node, block);
        insert(node, block, new_gain);
}

inline Gain tabu_bucket_queue::getKey(NodeID node, PartitionID block) {
        return m_gains->get_xy(node, block);
}
  
inline void tabu_bucket_queue::deleteNode(NodeID node, PartitionID block) {
        Count in_bucket_idx = m_queue_index->get_xy(node, block); 
        Gain  old_gain = m_gains->get_xy(node, block);
        unsigned address = old_gain + m_gain_span;

        if( m_buckets[address].size() > 1 ) {
                //swap current element with last element and pop_back
                std::pair< NodeID, PartitionID > p = m_buckets[address].back();

                m_queue_index->set_xy(p.first, p.second, in_bucket_idx) ; // update helper structure
                swap(m_buckets[address][in_bucket_idx], m_buckets[address].back());
                m_buckets[address].pop_back();
        } else {
                //size is 1
                m_buckets[address].pop_back();
                if( address == m_max_idx ) {
                        //update m_max_idx
                        while( m_max_idx != 0 )  {
                                m_max_idx--;
                                if(m_buckets[m_max_idx].size() > 0) {
                                        break;
                                }
                        }

                }
        }

        m_elements--;
        m_queue_index->set_xy(node, block, NOTINQUEUE);
        m_gains->set_xy(node, block, NOTINQUEUE);
}

inline bool tabu_bucket_queue::contains(NodeID node, PartitionID block) {
        return m_queue_index->get_xy(node, block) != NOTINQUEUE;
}


#endif /* end of include guard: BUCKET_PQ_EM8YJPA9 */
