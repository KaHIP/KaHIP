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
                void clear() { mapping_type.clear(); };
                NodeWeight & operator[](NodeID node) {return mapping_type[node];};

        private:
                linear_probing_hashmap mapping_type;
                PPartitionConfig m_config;
};

template <>
class hmap_wrapper <std::unordered_map<NodeID, NodeWeight> > {
        public:

                hmap_wrapper(PPartitionConfig & config) {
                        m_config = config;
                };

                virtual ~hmap_wrapper() {};

                void init( NodeID max_fill_count )  {};
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
                        // for large k, remember the entries that may be nonzero and reset only those
                        m_track = m_config.k > 32;
                        m_touched.clear();
                };

                void clear() { 
                        if( m_track && m_touched.size() < m_config.k ) {
                                for( ULONG i = 0; i < m_touched.size(); i++) {
                                        mapping_type[m_touched[i]] = 0;
                                }
                        } else {
                                for( ULONG k = 0; k < m_config.k; k++) {
                                        mapping_type[k] = 0;
                                }
                        }
                        m_touched.clear();
                };

                NodeWeight & operator[](NodeID node) {
                        if( m_track && mapping_type[node] == 0 ) m_touched.push_back(node);
                        return mapping_type[node];
                };

         private:
                std::vector<NodeWeight> mapping_type;
                std::vector<NodeID> m_touched;
                bool m_track;
                PPartitionConfig m_config;
};



#endif /* end of include guard: HMAP_WRAPPER_RQFK3ARC */
