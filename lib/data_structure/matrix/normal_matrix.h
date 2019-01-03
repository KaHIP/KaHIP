/******************************************************************************
 * normal_matrix.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef NORMAL_MATRIX_DAUJ4JMM
#define NORMAL_MATRIX_DAUJ4JMM

#include <vector>
#include <iostream>
#include "matrix.h"

class normal_matrix : public matrix {
public:
        normal_matrix(unsigned int dim_x, unsigned int dim_y, int lazy_init_val = 0) : m_dim_x (dim_x), 
                                                                                       m_dim_y (dim_y), 
                                                                                       m_lazy_init_val ( lazy_init_val ) {
                m_internal_matrix.resize(m_dim_x); //allocate the rest lazy
        };
        virtual ~normal_matrix() {};

        inline int get_xy(unsigned int x, unsigned int y) {
                if( m_internal_matrix[x].size() == 0 ) { 
                        return m_lazy_init_val;
                }
                return m_internal_matrix[x][y];
        };

        inline void set_xy(unsigned int x, unsigned int y, int value) {
                //resize the fields lazy
                if( m_internal_matrix[x].size() == 0 ) { 
                        m_internal_matrix[x].resize(m_dim_y);
                        for( unsigned y_1 = 0; y_1 < m_dim_y; y_1++) {
                                m_internal_matrix[x][y_1] = m_lazy_init_val;
                        }
                }
                m_internal_matrix[x][y] = value;
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
        std::vector< std::vector<int> > m_internal_matrix;
        unsigned int m_dim_x, m_dim_y;
        int m_lazy_init_val;
};


#endif /* end of include guard: NORMAL_MATRIX_DAUJ4JMM */
