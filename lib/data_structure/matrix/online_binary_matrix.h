/******************************************************************************
 * online_binary_matrix.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef ONLINE_BINARY_MATRIX_DAUJ4JMM
#define ONLINE_BINARY_MATRIX_DAUJ4JMM

#include <vector>
#include <iostream>
#include "matrix.h"

class online_binary_matrix : public matrix {
public:
        online_binary_matrix(unsigned int dim_x, unsigned int dim_y) : m_dim_x (dim_x), 
                                                                         m_dim_y (dim_y) {
        };
        
        void setPartitionConfig( PartitionConfig & config ) {
                int groups_size = config.group_sizes.size();
                config.compact_bin_id = new std::vector<unsigned int> (config.k,0);
                std::vector< int > interval_sizes(groups_size,0);
                interval_sizes[0] =  config.group_sizes[0]; 
                // for( unsigned i = 1; i < groups_size; i++) {
                //         interval_sizes[i] = config.group_sizes[i]*interval_sizes[i-1];
                // }
                // config.bit_sec_len = ceil(log2(config.k));

                config.bit_sec_len = 1;
                for( unsigned k = 0; k < groups_size; k++) {
                        int tmp = ceil(log2(config.group_sizes[k]));
                        if (tmp > config.bit_sec_len) {
                                config.bit_sec_len = tmp;
                        }
                }

                for (unsigned i = 0; i < config.k; i++) {
                        unsigned int lay_id = i;
                        for(int k=0; k < groups_size; k++) {
                                int remainder = lay_id % config.group_sizes[k];
                                lay_id = lay_id / config.group_sizes[k];
                                (*config.compact_bin_id)[i] += remainder << (k*config.bit_sec_len);
                        }
                }

                // for (unsigned i = 0; i < config.k; i++) {
                //         int k = groups_size-1;
                //         for(;k >= 0; k--) {
                //                 unsigned long long int lay_id = (unsigned long long int)i / interval_sizes[k];
                //                 (*config.compact_bin_id)[i] += lay_id << (k*config.bit_sec_len);
                //         }
                // }
 
                this->config = config;
        }

        virtual ~online_binary_matrix() {};

        inline int get_xy(unsigned int x, unsigned int y) {
                //now depending on x and y, generate distance
                int k = 0;
                unsigned long long int xor_x_y = (*config.compact_bin_id)[x] ^ (*config.compact_bin_id)[y];
                int count_leading_zeros = __builtin_clzll(xor_x_y);
                int total_n_bits = 8*sizeof(unsigned long long int);
                int clz = total_n_bits - count_leading_zeros -1;
                if (clz >= 0) {
                        k = (int)floor(clz / config.bit_sec_len);                      
                        return config.distances[k];
                } else  {
                        return 0;
                }       
        };

        inline void set_xy(unsigned int x, unsigned int y, int value) {
                // do nothing -- matrix cannot be modified
        };

        inline unsigned int get_x_dim() {return m_dim_x;};
        inline unsigned int get_y_dim() {return m_dim_y;};

        void print() {
                for( unsigned int i = 0; i < get_x_dim(); i++) {
                        for( unsigned int j = 0; j < get_y_dim(); j++) {
                                std::cout <<  get_xy(i,j) << " ";
                        }
                        std::cout <<  ""  << std::endl;
                }
        }

private:
        PartitionConfig config;
        unsigned int m_dim_x, m_dim_y;
};


#endif /* end of include guard: NORMAL_MATRIX_DAUJ4JMM */
