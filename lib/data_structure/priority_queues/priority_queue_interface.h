/******************************************************************************
 * priority_queue_interface.h 
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

#ifndef PRIORITY_QUEUE_INTERFACE_20ZSYG7R
#define PRIORITY_QUEUE_INTERFACE_20ZSYG7R

#include "definitions.h"

class priority_queue_interface {
        public:
                priority_queue_interface( ) {};
                virtual ~priority_queue_interface() {};

                /* returns the size of the priority queue */
                virtual NodeID size() = 0;
                virtual bool empty()  = 0 ;
               
                virtual void insert(NodeID id, Gain gain) = 0; 
                
                virtual Gain maxValue()     = 0;
                virtual NodeID maxElement() = 0;
                virtual NodeID deleteMax()  = 0;

                virtual void decreaseKey(NodeID node, Gain newGain) = 0;
                virtual void increaseKey(NodeID node, Gain newKey)  = 0;

                virtual void changeKey(NodeID element, Gain newKey) = 0;
                virtual Gain getKey(NodeID element)  = 0;
                virtual void deleteNode(NodeID node) = 0;
                virtual bool contains(NodeID node)   = 0;
};

typedef priority_queue_interface refinement_pq;

#endif /* end of include guard: PRIORITY_QUEUE_INTERFACE_20ZSYG7R */

