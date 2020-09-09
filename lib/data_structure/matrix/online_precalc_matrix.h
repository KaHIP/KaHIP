/******************************************************************************
 * online_precalc_matrix.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef ONLINE_PRECALC_MATRIX_DAUJ4JMM
#define ONLINE_PRECALC_MATRIX_DAUJ4JMM

#include <vector>
#include <iostream>
#include "matrix.h"

class online_precalc_matrix : public matrix {
public:
        online_precalc_matrix(unsigned int dim_x, unsigned int dim_y) : m_dim_x (dim_x), 
                                                                         m_dim_y (dim_y) {
        };
        
        void setPartitionConfig( PartitionConfig & config ) {
                int groups_size = config.group_sizes.size();
                config.bin_id = new std::vector<std::vector<int>> (config.k,std::vector<int>(groups_size));
                std::vector< int > interval_sizes(groups_size,0);
                interval_sizes[0] =  config.group_sizes[0]; 
                for( unsigned i = 1; i < groups_size; i++) {
                        interval_sizes[i] = config.group_sizes[i]*interval_sizes[i-1];
                }
                for (unsigned i = 0; i < config.k; i++) {
                        int k = groups_size-1;
                        for(;k >= 0; k--) {
                                (*config.bin_id)[i][k] = (unsigned)i / interval_sizes[k];
                        }
                }

                this->config = config;
        }

        virtual ~online_precalc_matrix() {};

        inline int get_xy(unsigned int x, unsigned int y) {
                //now depending on x and y, generate distance
                int k = config.group_sizes.size()-1;
                for(;k >= 0; k--) {
                        if( (*config.bin_id)[x][k] != (*config.bin_id)[y][k] ) {
                                break;
                        }
                }
                k++;
                return config.distances[k];
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
