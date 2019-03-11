/******************************************************************************
 * bucket_pq.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef BUCKET_PQ_EM8YJPA9
#define BUCKET_PQ_EM8YJPA9

#include <limits>
#include <unordered_map>
#include <utility>

#include "priority_queue_interface.h"

class bucket_pq : public priority_queue_interface {
        public:
                bucket_pq( const EdgeWeight & gain_span );

                ~bucket_pq() override = default;

                NodeID size() override;
                void insert(NodeID id, Gain gain) override;
                bool empty() override;

                Gain maxValue() override;
                NodeID maxElement() override;
                NodeID deleteMax() override;

                void decreaseKey(NodeID node, Gain newGain) override;
                void increaseKey(NodeID node, Gain newGain) override;

                void changeKey(NodeID element, Gain newKey) override;
                Gain getKey(NodeID element) override;
                void deleteNode(NodeID node) override;

                bool contains(NodeID node) override;

              private:
                NodeID     m_elements;
                EdgeWeight m_gain_span;
                unsigned   m_max_idx; //points to the non-empty bucket with the largest gain

                std::unordered_map<NodeID, std::pair<Count, Gain> > m_queue_index;
                std::vector< std::vector<NodeID> >             m_buckets;
};

inline bucket_pq::bucket_pq( const EdgeWeight & gain_span_input ) {
        m_elements  = 0;
        m_gain_span = gain_span_input;
        m_max_idx   = 0;

        m_buckets.resize(2*m_gain_span+1);
}

inline NodeID bucket_pq::size() {
        return m_elements;
}

inline void bucket_pq::insert(NodeID node, Gain gain) {
        unsigned address = gain + m_gain_span;
        if(address > m_max_idx) {
                m_max_idx = address;
        }

        m_buckets[address].push_back( node );
        m_queue_index[node].first  = m_buckets[address].size() - 1; //store position
        m_queue_index[node].second = gain;

        m_elements++;
}

inline bool bucket_pq::empty( ) {
        return m_elements == 0;
}

inline Gain bucket_pq::maxValue( ) {
        return m_max_idx - m_gain_span;
}

inline NodeID bucket_pq::maxElement( ) {
        return m_buckets[m_max_idx].back();
}

inline NodeID bucket_pq::deleteMax() {
       NodeID node = m_buckets[m_max_idx].back();
       m_buckets[m_max_idx].pop_back();
       m_queue_index.erase(node);

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
       return node;
}

inline void bucket_pq::decreaseKey(NodeID node, Gain new_gain) {
        changeKey( node, new_gain );
}

inline void bucket_pq::increaseKey(NodeID node, Gain new_gain) {
        changeKey( node, new_gain );
}

inline Gain bucket_pq::getKey(NodeID node) {
        return m_queue_index[node].second;
}

inline void bucket_pq::changeKey(NodeID node, Gain new_gain) {
        deleteNode(node);
        insert(node, new_gain);
}

inline void bucket_pq::deleteNode(NodeID node) {
        ASSERT_TRUE(m_queue_index.find(node) != m_queue_index.end());
        Count in_bucket_idx = m_queue_index[node].first;
        Gain  old_gain      = m_queue_index[node].second;
        unsigned address    = old_gain + m_gain_span;

        if( m_buckets[address].size() > 1 ) {
                //swap current element with last element and pop_back
                m_queue_index[m_buckets[address].back()].first = in_bucket_idx; // update helper structure
                std::swap(m_buckets[address][in_bucket_idx], m_buckets[address].back());
                m_buckets[address].pop_back();
        } else {
                //size is 1
                m_buckets[address].pop_back();
                if( address == m_max_idx ) {
                        //update max_idx
                        while( m_max_idx != 0 )  {
                                m_max_idx--;
                                if(m_buckets[m_max_idx].size() > 0) {
                                        break;
                                }
                        }

                }
        }

        m_elements--;
        m_queue_index.erase(node);
}

inline bool bucket_pq::contains(NodeID node) {
        return m_queue_index.find(node) != m_queue_index.end();
}


#endif /* end of include guard: BUCKET_PQ_EM8YJPA9 */
