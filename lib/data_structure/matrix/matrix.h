/******************************************************************************
 * matrix.h 
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

#ifndef MATRIX_BHHZ9T7P
#define MATRIX_BHHZ9T7P

class matrix {
public:
        matrix(unsigned int dim_x, unsigned int dim_y) {};
        matrix() {};
        virtual ~matrix() {};

        virtual int  get_xy(unsigned int x, unsigned int y)            = 0;
        virtual void set_xy(unsigned int x, unsigned int y, int value) = 0;
};


#endif /* end of include guard: MATRIX_BHHZ9T7P */
