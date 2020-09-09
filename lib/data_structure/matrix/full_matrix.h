/******************************************************************************
 * full_matrix.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef FULL_MATRIX_DAUJ4JMM
#define FULL_MATRIX_DAUJ4JMM

#include <vector>
#include <iostream>
#include "matrix.h"

#define MM(a,b) m_internal_matrix[((a)*m_dim_x)+(b)]

class full_matrix : public matrix {
public:
        full_matrix(unsigned int dim_x, unsigned int dim_y, int lazy_init_val = 0) : m_dim_x (dim_x), 
                                                                                       m_dim_y (dim_y), 
                                                                                       m_lazy_init_val ( lazy_init_val ) {
                m_internal_matrix.resize(m_dim_x*m_dim_x); 
        };
        virtual ~full_matrix() {};

        void setPartitionConfig( PartitionConfig & config ) {
                std::vector< int > interval_sizes;
                // this->config = config;
                interval_sizes.resize(config.group_sizes.size());
                interval_sizes[0] = config.group_sizes[0]; 
                for( unsigned i = 1; i < interval_sizes.size(); i++) {
                        interval_sizes[i] = config.group_sizes[i]*interval_sizes[i-1];
                }

                for (int x = 0; x< m_dim_x; x++ ) {
                        for( unsigned y = 0; y < m_dim_y; y++) {
                                if (x == y) {
                                        MM(x,y) = 0;
                                        continue;
                                }

                                int k = config.group_sizes.size()-1;
                
                                for(;k >= 0; k--) {
                                        int interval_a = x / interval_sizes[k];
                                        int interval_b = y / interval_sizes[k];
                                        if( interval_a != interval_b ) {
                                                break;
                                        }
                                }
                                k++;
                                MM(x,y) = config.distances[k];
                        }
                }
        }

        inline int get_xy(unsigned int x, unsigned int y) {
                // int k = config.group_sizes.size()-1;

                // for(;k >= 0; k--) {
                //         int interval_a = x / interval_sizes[k];
                //         int interval_b = y / interval_sizes[k];
                //         if( interval_a != interval_b ) {
                //                 break;
                //         }
                // }
                // k++;

                // if (MM(x,y) != config.distances[k]) {
                //         std::cout << "Difference in Result ! (MARCELO)\n" ;
                // }

                return MM(x,y);
        };

        inline void set_xy(unsigned int x, unsigned int y, int value) {

        }

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
        // PartitionConfig config;

        std::vector<int> m_internal_matrix;
        unsigned int m_dim_x, m_dim_y;
        int m_lazy_init_val;
        // std::vector< int > interval_sizes;
};


#endif /* end of include guard: FULL_MATRIX_DAUJ4JMM */
