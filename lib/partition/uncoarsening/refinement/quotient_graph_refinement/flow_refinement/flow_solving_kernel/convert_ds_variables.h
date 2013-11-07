/******************************************************************************
 * convert_ds_variables.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * based on HiPr  -- Copyright 1995, 2000 by IG Systems, Inc., igsys@eclipse.net 
 * "Maximum flow - highest lavel push-relabel algorithm"
 * Copyright for this release under gpl 3.0 granted by Andrew Goldberg
 *
 * Comment: we used the implementation of hipr and put them into a header file here
 *
 ******************************************************************************
 * Copyright (C) 2013 Christian Schulz <christian.schulz@kit.edu>
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

	long    n,                      /* internal number of nodes */
        node_min = 0,                     /* minimal no of node  */
        node_max = 0,                     /* maximal no of nodes */
        *arc_first = NULL,                /* internal array for holding
                                             - node degree
                                             - position of the first outgoing arc */
        *arc_tail = NULL,                  /* internal array: tails of the arcs */
        source = 0,                        /* no of the source */
        sink = 0,                          /* no of the sink */
        /* temporary variables carrying no of nodes */
        head, tail, i;

        long    m,                       /* internal number of arcs */
                /* temporary variables carrying no of arcs */
                last, arc_num, arc_new_num;

        node    *nodes = NULL,             /* pointer to the node structure */
                *head_p,
                *ndp;

        arc     *arcs = NULL,              /* pointer to the arc structure */
                *arc_current = NULL,
                *arc_new,
                *arc_tmp;

        long    *acap = NULL,              /* array of capasities */
                cap;                     /* capasity of the current arc */

        long pos_current = 0;               /* [> 2*no_alines <] */


