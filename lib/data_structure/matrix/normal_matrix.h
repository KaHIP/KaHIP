/******************************************************************************
 * normal_matrix.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#ifndef NORMAL_MATRIX_DAUJ4JMM
#define NORMAL_MATRIX_DAUJ4JMM

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

private:
        std::vector< std::vector<int> > m_internal_matrix;
        unsigned int m_dim_x, m_dim_y;
        int m_lazy_init_val;
};


#endif /* end of include guard: NORMAL_MATRIX_DAUJ4JMM */
