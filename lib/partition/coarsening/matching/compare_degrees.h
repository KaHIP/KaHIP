/******************************************************************************
 * compare_degrees.h 
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

#ifndef COMPARE_DEGREES_750FUZ7Z
#define COMPARE_DEGREES_750FUZ7Z

class compare_degrees : public std::binary_function<EdgeWeight, EdgeWeight, bool> {
        public:
                compare_degrees(std::vector<EdgeWeight> * degrees) : m_node_degrees(degrees) {};
                virtual ~compare_degrees() {};

                bool operator() (const EdgeWeight left, const EdgeWeight right ) {
                        return (*m_node_degrees)[left] < (*m_node_degrees)[right];
                }

        private:
                std::vector<EdgeWeight> * m_node_degrees;
};


#endif /* end of include guard: COMPARE_DEGREES_750FUZ7Z */
