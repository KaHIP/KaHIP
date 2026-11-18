/******************************************************************************
 * hmap_wrapper.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef HMAP_WRAPPER_RQFK3ARC
#define HMAP_WRAPPER_RQFK3ARC

#include "data_structure/linear_probing_hashmap.h"

template <typename T>
class hmap_wrapper {
        public:

                hmap_wrapper(PPartitionConfig & config) {
                        m_config = config;
                };

                virtual ~hmap_wrapper() {};

                void init( NodeID max_fill_count );
                void prepare( NodeID degree );
                void clear();
                NodeWeight & operator[](NodeID node);

private:
                T mapping_type;
                PPartitionConfig m_config;
};

template <>
class hmap_wrapper < linear_probing_hashmap > {
        public:

                hmap_wrapper(PPartitionConfig & config) {
                        m_config = config;
                };

                virtual ~hmap_wrapper() {};

                void init(NodeID max_fill_count )  {mapping_type.init(max_fill_count, m_config.ht_fill_factor);};

                // a node of small degree keeps its neighbor labels in a short array that is searched
                // linearly instead of the hash table (the weights per label are the same)
                void prepare( NodeID degree ) { m_small = degree <= SMALL_DEGREE; m_small_count = 0; };

                void clear() { 
                        if( m_small ) m_small_count = 0;
                        else          mapping_type.clear(); 
                };

                NodeWeight & operator[](NodeID node) {
                        if( m_small ) {
                                for( unsigned i = 0; i < m_small_count; i++) {
                                        if( m_small_keys[i] == node ) return m_small_values[i];
                                }
                                m_small_keys[m_small_count]   = node;
                                m_small_values[m_small_count] = 0;
                                return m_small_values[m_small_count++];
                        }
                        return mapping_type[node];
                };

        private:
                static const unsigned SMALL_DEGREE = 16;
                linear_probing_hashmap mapping_type;
                PPartitionConfig m_config;
                bool m_small = false;
                unsigned m_small_count = 0;
                NodeID m_small_keys[SMALL_DEGREE];
                NodeWeight m_small_values[SMALL_DEGREE];
};

template <>
class hmap_wrapper <std::unordered_map<NodeID, NodeWeight> > {
        public:

                hmap_wrapper(PPartitionConfig & config) {
                        m_config = config;
                };

                virtual ~hmap_wrapper() {};

                void init( NodeID max_fill_count )  {};
                void prepare( NodeID degree ) {};
                void clear() { mapping_type.clear(); };
                NodeWeight & operator[](NodeID node) {return mapping_type[node];};

        private:
                std::unordered_map<NodeID, NodeWeight> mapping_type;
                PPartitionConfig m_config;
};

template <>
class hmap_wrapper <std::vector<NodeWeight> > {
        public:

                hmap_wrapper(PPartitionConfig & config) {
                        m_config = config;
                };

                virtual ~hmap_wrapper() {};

                void init( NodeID max_fill )  {
                        mapping_type.resize(m_config.k);
                        for( ULONG k = 0; k < m_config.k; k++) {
                                mapping_type[k] = 0;
                        }

                };

                void prepare( NodeID degree ) {};

                // reset only the entries used since the last clear (an entry that is still 0 is
                // recorded again on its next use, which is harmless)
                void clear() { 
                        for( size_t i = 0; i < m_touched.size(); i++) {
                                mapping_type[m_touched[i]] = 0;
                        }
                        m_touched.clear();
                };

                NodeWeight & operator[](NodeID node) {
                        NodeWeight & value = mapping_type[node];
                        if( value == 0 ) m_touched.push_back(node);
                        return value;
                };

         private:
                std::vector<NodeWeight> mapping_type;
                std::vector<NodeID> m_touched;
                PPartitionConfig m_config;
};



#endif /* end of include guard: HMAP_WRAPPER_RQFK3ARC */
