/******************************************************************************
 * only_assign.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
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
 
