/******************************************************************************
 * construct_distance_matrix.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Graph Partitioning 
 ******************************************************************************
 * Copyright (C) 2017 Christian Schulz 
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

#ifndef CONSTRUCT_DISTANCE_MATRIX_HUPBUT8O
#define CONSTRUCT_DISTANCE_MATRIX_HUPBUT8O

#include "data_structure/matrix/matrix.h"
#include "partition_config.h"
#include "tools/random_functions.h"

class construct_distance_matrix {
public:
        construct_distance_matrix();
        virtual ~construct_distance_matrix();

        void construct_matrix( PartitionConfig & config, matrix & D ) {
                //check wether distance matrix is a square matrix
                if(D.get_x_dim() != D.get_y_dim()) {
                        std::cout <<  "distance matrix is not symmetric."  << std::endl;
                        exit(0);
                }

                switch( config.distance_construction_algorithm ) {
                        case DIST_CONST_RANDOM:
                                construct_matrix_random( config, D);
                                break;
                        case DIST_CONST_IDENTITY:
                                construct_matrix_identity( config, D);
                                break;
                        case DIST_CONST_HIERARCHY: 
                                construct_matrix_hierarchy( config, D);
                                break;
                        case DIST_CONST_HIERARCHY_ONLINE: 
                                break;
                        default: 
                                construct_matrix_random( config, D );
                }
        };

private:

        void construct_matrix_random( PartitionConfig & config, matrix & D ) {
                for( unsigned int i = 0; i < D.get_x_dim(); i++) {
                        for( unsigned int j = 0; j <= i ; j++) {
                                NodeWeight value = random_functions::nextInt(1,100);
                                D.set_xy(i,j, value);
                                D.set_xy(j,i, value);
                        }
                }
        }

        void construct_matrix_identity( PartitionConfig & config, matrix & D ) {
                for( unsigned int i = 0; i < D.get_x_dim(); i++) {
                        for( unsigned int j = 0; j <= i ; j++) {
                                D.set_xy(i,j, 1);
                                D.set_xy(j,i, 1);
                        }
                }
        }

        void construct_matrix_hierarchy( PartitionConfig & config, matrix & D ) {
                std::vector< int > interval_sizes(config.group_sizes.size(),0);
                interval_sizes[0] =  config.group_sizes[0]; 
                for( unsigned i = 1; i < interval_sizes.size(); i++) {
                        interval_sizes[i] = config.group_sizes[i]*interval_sizes[i-1];
                }
                
                PRINT(std::cout <<  "total num cores " << interval_sizes[interval_sizes.size()-1]  << std::endl;)
                for( unsigned int i = 0; i < D.get_x_dim(); i++) {
                        for( unsigned int j = 0; j <= i ; j++) {
                                //now depending on i and j, generate distance
                                int k = config.group_sizes.size()-1;
                                for(;k >= 0; k--) {
                                        int interval_a = i / interval_sizes[k];
                                        int interval_b = j / interval_sizes[k];
                                        if( interval_a != interval_b ) {
                                                break;
                                        }
                                }
                                k++;

                                NodeWeight distance = config.distances[k];
                                D.set_xy(i,j, distance);
                                D.set_xy(j,i, distance);
                        }
                }
        }
};


#endif /* end of include guard: CONSTRUCT_DISTANCE_MATRIX_HUPBUT8O */
