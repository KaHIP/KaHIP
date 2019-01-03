/******************************************************************************
 * convert_ds_variables.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
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


