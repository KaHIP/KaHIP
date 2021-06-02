/******************************************************************************
 * hmap_wrapper.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef HMAP_WRAPPER_RQFK3ARC
#define HMAP_WRAPPER_RQFK3ARC

#include "data_structure/linear_probing_hashmap.h"
#include "partition_config.h"

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
                hmap_wrapper() {
                };
                hmap_wrapper(PPartitionConfig & config) {
                        m_config = config;
                };

                virtual ~hmap_wrapper() {};

                void init(NodeID max_fill_count )  {
                        mapping_type.clear();
                        mapping_type.reserve(max_fill_count);
                        mapping_type_small.assign(max_fill_count, 0);
                        small_thresold = max_fill_count;
                };
                void clear() { 
                        mapping_type.clear(); 
                        while (mapping_type_small.size() > small_thresold) {
                                mapping_type_small[mapping_type_small.back()] = 0;
                                mapping_type_small.pop_back();
                        }
                };
                NodeWeight & operator[](NodeID node) {
                        if (node < small_thresold) {
                                if (!mapping_type_small[node]) {
                                        mapping_type_small.emplace_back(node);
                                }
                                return mapping_type_small[node];
                        }
                        return mapping_type[node];
                };
                void erase(NodeID node) {
                        if (node < small_thresold) {
                                mapping_type_small[node] = 0;
                                return;
                        }
                        mapping_type.erase(node);
                }

        private:
                size_t small_thresold;
                std::vector<NodeWeight> mapping_type_small;
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

                void clear() { 
                        for( ULONG k = 0; k < m_config.k; k++) {
                                mapping_type[k] = 0;
                        }
                };

                NodeWeight & operator[](NodeID node) {return mapping_type[node];};

         private:
                std::vector<NodeWeight> mapping_type;
                PPartitionConfig m_config;
};



#endif /* end of include guard: HMAP_WRAPPER_RQFK3ARC */
