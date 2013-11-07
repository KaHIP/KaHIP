/******************************************************************************
 * linear_ordering_n_assign.h 
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


        /********** ordering arcs - linear time algorithm ***********/

        /* first arc from the first node */
        ( nodes + node_min ) -> first = arcs;

        /* before below loop arc_first[i+1] is the number of arcs outgoing from i;
           after this loop arc_first[i] is the position of the first 
           outgoing from node i arcs after they would be ordered;
           this value is transformed to pointer and written to node.first[i]
           */

        for ( i = node_min + 1; i <= node_max + 1; i ++ ) 
        {
                arc_first[i]          += arc_first[i-1];
                ( nodes + i ) -> first = arcs + arc_first[i];
        }


        for ( i = node_min; i < node_max; i ++ ) /* scanning all the nodes  
                                                    exept the last*/
        {

                last = ( ( nodes + i + 1 ) -> first ) - arcs;
                /* arcs outgoing from i must be cited    
                   from position arc_first[i] to the position
                   equal to initial value of arc_first[i+1]-1  */

                for ( arc_num = arc_first[i]; arc_num < last; arc_num ++ )
                { tail = arc_tail[arc_num];

                        while ( tail != i )
                                /* the arc no  arc_num  is not in place because arc cited here
                                   must go out from i;
                                   we'll put it to its place and continue this process
                                   until an arc in this position would go out from i */

                        { arc_new_num  = arc_first[tail];
                                arc_current  = arcs + arc_num;
                                arc_new      = arcs + arc_new_num;

                                /* arc_current must be cited in the position arc_new    
                                   swapping these arcs:                                 */

                                head_p               = arc_new -> head;
                                arc_new -> head      = arc_current -> head;
                                arc_current -> head  = head_p;

                                cap                 = arc_new -> resCap;
                                arc_new -> resCap     = arc_current -> resCap;
                                arc_current -> resCap = cap;

                                if ( arc_new != arc_current -> rev )
                                {
                                        arc_tmp                = arc_new -> rev;
                                        arc_new  -> rev     = arc_current -> rev;
                                        arc_current -> rev  = arc_tmp;

                                        ( arc_current -> rev ) -> rev = arc_current;
                                        ( arc_new     -> rev ) -> rev = arc_new;
                                }

                                arc_tail[arc_num] = arc_tail[arc_new_num];
                                arc_tail[arc_new_num] = tail;

                                /* we increase arc_first[tail]  */
                                arc_first[tail] ++ ;

                                tail = arc_tail[arc_num];
                        }
                }
                /* all arcs outgoing from  i  are in place */
        }       

        /* -----------------------  arcs are ordered  ------------------------- */
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
  
