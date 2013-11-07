/******************************************************************************
 * only_assign.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * based on HiPr  -- Copyright 1995, 2000 by IG Systems, Inc., igsys@eclipse.net 
 * "Maximum flow - highest lavel push-relabel algorithm"
 * Copyright for this release under gpl 3.0 granted by Andrew Goldberg
 *
 * Comment: we used the implementation of hipr and put them into a header file here
 *
 ******************************************************************************
 * Copyright (modifications) (C) 2013 Christian Schulz <christian.schulz@kit.edu>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
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

        /*----------- constructing lists ---------------*/


        for ( ndp = nodes + node_min; ndp <= nodes + node_max;  ndp ++ )
                ndp -> first = (arc*) NULL;

        for ( arc_current = arcs + (2*m-1); arc_current >= arcs; arc_current -- )
        {
                arc_num = arc_current - arcs;
                tail = arc_tail [arc_num];
                ndp = nodes + tail;
                /* avg
                   arc_current -> next = ndp -> first;
                   */
                ndp -> first = arc_current;
        }


        /* ----------- assigning output values ------------*/
        *m_ad = m;
        *n_ad = node_max - node_min + 1;
        *source_ad = nodes + source;
        *sink_ad   = nodes + sink;
        *node_min_ad = node_min;
        *nodes_ad = nodes + node_min;
        *arcs_ad = arcs;
        *cap_ad = acap;

        for ( arc_current = arcs, arc_num = 0; 
                        arc_num < 2*m;
                        arc_current ++, arc_num ++
            )
                acap [ arc_num ] = arc_current -> resCap; 

        /* free internal memory */
        free ( arc_first ); free ( arc_tail );
        free_nodes = nodes;
 
