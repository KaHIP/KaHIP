/******************************************************************************
 * flow_solver.cpp 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * based on HiPr  -- Copyright 1995, 2000 by IG Systems, Inc., igsys@eclipse.net 
 * "Maximum flow - highest lavel push-relabel algorithm"
 * Copyright for this release gpl 3.0 granted by Andrew Goldberg
 *
 * comment: we used the implementation of hipr and put them into a class structure here
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



#define CHECK_SOLUTION
#define PRINT_CUT
#define CUT_ONLY 

#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <sstream>

#include "flow_solver.h"
#include "flow_macros.h"
#include "most_balanced_minimum_cuts/most_balanced_minimum_cuts.h"


flow_solver::flow_solver() {
        pushCnt         = 0;       /* number of pushes */
        relabelCnt      = 0;       /* number of internal_relabels */
        updateCnt       = 0;       /* number of updates */
        gapCnt          = 0;       /* number of internal_gaps */
        gNodeCnt        = 0;       /* number of nodes after internal_gap */
        workSinceUpdate = 0;       /* the number of arc scans since last update */
        nodes           = NULL;
        arcs            = NULL;
        cap             = NULL;
        buckets         = NULL;
        free_nodes      = NULL;
}

flow_solver::~flow_solver() {
        free(arcs);
        free(cap);
        free(buckets);
        free(free_nodes);
}


void flow_solver::internal_stage_one() {

        node   *i;
        bucket  *l;             /* current bucket */


#if defined(INIT_UPDATE) || defined(OLD_INIT) || defined(WAVE_INIT)
        internal_global_update ();
#endif

        workSinceUpdate = 0;

#ifdef WAVE_INIT
        internal_wave();
#endif  

        /* main loop */
        while ( aMax >= aMin ) {
                l = buckets + aMax;
                i = l->firstActive;

                if (i == sentinelNode) {
                        aMax--;
                }
                else {
                        aRemove(l,i);
                        assert(i->excess > 0);
                        internal_discharge (i);

                        if (aMax < aMin)
                                break;

                        /* is it time for global update? */
                        if (workSinceUpdate * globUpdtFreq > nm) {
                                internal_global_update ();
                                workSinceUpdate = 0;
                        }

                }

        } /* end of the main loop */

        flow = sink -> excess;
} 


void flow_solver::internal_stage_two() {
        node *i, *j, *tos, *bos, *restart, *r;
        arc *a;
        cType delta;

        /* deal with self-loops */
        forAllNodes(i) {
                forAllArcs(i,a)
                        if ( a -> head == i ) {
                                a -> resCap = cap[a - arcs];
                        }
        }

        /* initialize */
        tos = bos = NULL;
        forAllNodes(i) {
                i -> d = WHITE;
                //    buckets[i-nodes].firstActive = NULL;
                buckets[i-nodes].firstActive = sentinelNode;
                i -> current = i -> first;
        }

        /* eliminate flow cycles, topologicaly order vertices */
        forAllNodes(i)
                if (( i -> d == WHITE ) && ( i -> excess > 0 ) &&
                                ( i != source ) && ( i != sink )) {
                        r = i;
                        r -> d = GREY;
                        do {
                                for ( ; i->current != (i+1)->first; i->current++) {
                                        a = i -> current;
                                        if (( cap[a - arcs] == 0 ) && ( a -> resCap > 0 )) { 
                                                j = a -> head;
                                                if ( j -> d == WHITE ) {
                                                        /* start scanning j */
                                                        j -> d = GREY;
                                                        buckets[j-nodes].firstActive = i;
                                                        i = j;
                                                        break;
                                                }
                                                else
                                                        if ( j -> d == GREY ) {
                                                                /* find minimum flow on the cycle */
                                                                delta = a -> resCap;
                                                                while ( 1 ) {
                                                                        delta = min ( delta, j -> current -> resCap );
                                                                        if ( j == i )
                                                                                break;
                                                                        else
                                                                                j = j -> current -> head;
                                                                }

                                                                /* remove delta flow units */
                                                                j = i;
                                                                while ( 1 ) {
                                                                        a = j -> current;
                                                                        a -> resCap -= delta;
                                                                        a -> rev -> resCap += delta;
                                                                        j = a -> head;
                                                                        if ( j == i )
                                                                                break;
                                                                }

                                                                /* backup DFS to the first saturated arc */
                                                                restart = i;
                                                                for ( j = i -> current -> head; j != i; j = a -> head ) {
                                                                        a = j -> current;
                                                                        if (( j -> d == WHITE ) || ( a -> resCap == 0 )) {
                                                                                j -> current -> head -> d = WHITE;
                                                                                if ( j -> d != WHITE )
                                                                                        restart = j;
                                                                        }
                                                                }

                                                                if ( restart != i ) {
                                                                        i = restart;
                                                                        i->current++;
                                                                        break;
                                                                }
                                                        }
                                        }
                                }

                                if (i->current == (i+1)->first) {
                                        /* scan of i complete */
                                        i -> d = BLACK;
                                        if ( i != source ) {
                                                if ( bos == NULL ) {
                                                        bos = i;
                                                        tos = i;
                                                }
                                                else {
                                                        i -> bNext = tos;
                                                        tos = i;
                                                }
                                        }

                                        if ( i != r ) {
                                                i = buckets[i-nodes].firstActive;
                                                i->current++;
                                        }
                                        else
                                                break;
                                }
                        } while ( 1 );
                }


        /* return excesses */
        /* note that sink is not on the stack */
        if ( bos != NULL ) {
                for ( i = tos; i != bos; i = i -> bNext ) {
                        a = i -> first;
                        while ( i -> excess > 0 ) {
                                if (( cap[a - arcs] == 0 ) && ( a -> resCap > 0 )) {
                                        if (a->resCap < i->excess)
                                                delta = a->resCap;
                                        else
                                                delta = i->excess;
                                        a -> resCap -= delta;
                                        a -> rev -> resCap += delta;
                                        i -> excess -= delta;
                                        a -> head -> excess += delta;
                                }
                                a++;
                        }
                }
                /* now do the bottom */
                i = bos;
                a = i -> first;
                while ( i -> excess > 0 ) {
                        if (( cap[a - arcs] == 0 ) && ( a -> resCap > 0 )) {
                                if (a->resCap < i->excess)
                                        delta = a->resCap;
                                else
                                        delta = i->excess;
                                a -> resCap -= delta;
                                a -> rev -> resCap += delta;
                                i -> excess -= delta;
                                a -> head -> excess += delta;
                        }
                        a++;
                }
        }

}

void flow_solver::internal_global_update() {

        node  *i, *j;       /* node pointers */
        arc   *a;           /* current arc pointers  */
        bucket *l, *jL;          /* bucket */
        long curDist, jD;
        long state;


        updateCnt ++;

        /* initialization */

        forAllNodes(i)
                i -> d = n;
        sink -> d = 0;

        for (l = buckets; l <= buckets + dMax; l++) {
                l -> firstActive   = sentinelNode;
                l -> firstInactive  = sentinelNode;
        }

        dMax = aMax = 0;
        aMin = n;

        /* breadth first search */

        // add sink to bucket zero

        iAdd(buckets, sink);
        for (curDist = 0; 1; curDist++) {

                state = 0;
                l = buckets + curDist;
                jD = curDist + 1;
                jL = l + 1;
                /*
                   jL -> firstActive   = sentinelNode;
                   jL -> firstInactive  = sentinelNode;
                   */

                if ((l->firstActive == sentinelNode) && 
                        (l->firstInactive == sentinelNode))
                        break;

                while (1) {

                        switch (state) {
                                case 0: 
                                        i = l->firstInactive;
                                        state = 1;
                                        break;
                                case 1:
                                        i = i->bNext;
                                        break;
                                case 2:
                                        i = l->firstActive;
                                        state = 3;
                                        break;
                                case 3:
                                        i = i->bNext;
                                        break;
                                default: 
                                        assert(0);
                                        break;
                        }

                        if (i == sentinelNode) {
                                if (state == 1) {
                                        state = 2;
                                        continue;
                                }
                                else {
                                        assert(state == 3);
                                        break;
                                }
                        }

                        /* scanning arcs incident to node i */
                        forAllArcs(i,a) {
                                if (a->rev->resCap > 0 ) {
                                        j = a->head;
                                        if (j->d == n) {
                                                j->d = jD;
                                                j->current = j->first;
                                                if (jD > dMax) dMax = jD;

                                                if (j->excess > 0) {
                                                        /* put into active list */
                                                        aAdd(jL,j);
                                                }
                                                else {
                                                        /* put into inactive list */
                                                        iAdd(jL,j);
                                                }
                                        }
                                }
                        } /* node i is scanned */ 
                }
        }

} /* end of global update */


void flow_solver::internal_check_max()
{
        bucket *l;

        for (l = buckets + dMax + 1; l < buckets + n; l++) {
                assert(l->firstActive == sentinelNode);
                assert(l->firstInactive == sentinelNode);
        }
}

void flow_solver::internal_init( )
{
        node  *i;        /* current node */
        int overflowDetected;
        bucket *l;
        arc *a;
#ifdef EXCESS_TYPE_LONG
        double testExcess;
#endif
#ifndef OLD_INIT
        unsigned long delta;
#endif

        // initialize excesses

        forAllNodes(i) {
                i->excess = 0;
                i->current = i->first;
                forAllArcs(i, a)
                        a->resCap = cap[a-arcs];
        }

        for (l = buckets; l <= buckets + n-1; l++) {
                l -> firstActive   = sentinelNode;
                l -> firstInactive  = sentinelNode;
        }

        overflowDetected = 0;
#ifdef EXCESS_TYPE_LONG
        testExcess = 0;
        forAllArcs(source,a) {
                if (a->head != source) {
                        testExcess += a->resCap;
                }
        }
        if (testExcess > MAXLONG) {
                printf("c WARNING: excess overflow. See README for details.\nc\n");
                overflowDetected = 1;
        }
#endif
#ifdef OLD_INIT
        source -> excess = MAXLONG;
#else
        if (overflowDetected) {
                source -> excess = MAXLONG;
        }
        else {
                source->excess = 0;
                forAllArcs(source,a) {
                        if (a->head != source) {
                                pushCnt ++;
                                delta = a -> resCap;
                                a -> resCap -= delta;
                                (a -> rev) -> resCap += delta;
                                a->head->excess += delta;
                        }
                }
        }

        /*  setup labels and buckets */
        l = buckets + 1;

        aMax = 0;
        aMin = n;

        forAllNodes(i) {
                if (i == sink) {
                        i->d = 0;
                        iAdd(buckets,i);
                        continue;
                }
                if ((i == source) && (!overflowDetected)) {
                        i->d = n;
                }
                else
                        i->d = 1;
                if (i->excess > 0) {
                        /* put into active list */
                        aAdd(l,i);
                }
                else { /* i -> excess == 0 */
                        /* put into inactive list */
                        if (i->d < n)
                                iAdd(l,i);
                }
        }
        dMax = 1;
#endif

} /* end of init */

int flow_solver::internal_allocDS( )
{

        nm = ALPHA * n + m;
        /*
           queue = (node**) calloc ( n, sizeof (node*) );
           if ( queue == NULL ) return ( 1 );
           qLast = queue + n - 1;
           qInit();
           */
        buckets = (bucket*) calloc ( n+2, sizeof (bucket) );
        if ( buckets == NULL ) return ( 1 );

        sentinelNode = nodes + n;
        sentinelNode->first = arcs + 2*m;

        return ( 0 );

} /* end of allocate */

/* internal_gap internal_relabeling */

int flow_solver::internal_gap (bucket* emptyB)
{

        bucket *l;
        node  *i; 
        long  r;           /* index of the bucket before l  */
        int   cc;          /* cc = 1 if no nodes with positive excess before
                              the internal_gap */

        gapCnt ++;
        r = ( emptyB - buckets ) - 1;

        /* set labels of nodes beyond the internal_gap to "infinity" */
        for ( l = emptyB + 1; l <= buckets + dMax; l ++ ) {
                /* this does nothing for high level selection 
                   for (i = l -> firstActive; i != sentinelNode; i = i -> bNext) {
                   i -> d = n;
                   gNodeCnt++;
                   }
                   l -> firstActive = sentinelNode;
                   */

                for ( i = l -> firstInactive; i != sentinelNode; i = i -> bNext ) {
                        i -> d = n;
                        gNodeCnt ++;
                }

                l -> firstInactive = sentinelNode;
        }

        cc = ( aMin > r ) ? 1 : 0;

        dMax = r;
        aMax = r;

        return ( cc );

}

/*--- internal_relabelling node i */

long flow_solver::internal_relabel (node *i)
{

        node  *j;
        long  minD;     /* minimum d of a node reachable from i */
        arc   *minA;    /* an arc which leads to the node with minimal d */
        arc   *a;

        assert(i->excess > 0);

        relabelCnt++;
        workSinceUpdate += BETA;

        i->d = minD = n;
        minA = NULL;

        /* find the minimum */
        forAllArcs(i,a) {
                workSinceUpdate++;
                if (a -> resCap > 0) {
                        j = a -> head;
                        if (j->d < minD) {
                                minD = j->d;
                                minA = a;
                        }
                }
        }

        minD++;

        if (minD < n) {

                i->d = minD;
                i->current = minA;

                if (dMax < minD) dMax = minD;

        } /* end of minD < n */

        return ( minD );

} /* end of internal_relabel */


/* internal_discharge: push flow out of i until i becomes inactive */

void flow_solver::internal_discharge (node* i)
{

        node  *j;                 /* sucsessor of i */
        long  jD;                 /* d of the next bucket */
        bucket *lj;               /* j's bucket */
        bucket *l;                /* i's bucket */
        arc   *a;                 /* current arc (i,j) */
        cType  delta;
        arc *stopA;

        assert(i->excess > 0);
        assert(i != sink);
        do {

                jD = i->d - 1;
                l = buckets + i->d;

                /* scanning arcs outgoing from  i  */
                for (a = i->current, stopA = (i+1)->first; a != stopA; a++) {
                        if (a -> resCap > 0) {
                                j = a -> head;

                                if (j->d == jD) {
                                        pushCnt ++;
                                        if (a->resCap < i->excess)
                                                delta = a->resCap;
                                        else
                                                delta = i->excess;
                                        a->resCap -= delta;
                                        a->rev->resCap += delta;

                                        if (j != sink) {

                                                lj = buckets + jD;

                                                if (j->excess == 0) {
                                                        /* remove j from inactive list */
                                                        iDelete(lj,j);
                                                        /* add j to active list */
                                                        aAdd(lj,j);
                                                }
                                        }

                                        j -> excess += delta;
                                        i -> excess -= delta;

                                        if (i->excess == 0) break;

                                } /* j belongs to the next bucket */
                        } /* a  is not saturated */
                } /* end of scanning arcs from  i */

                if (a == stopA) {
                        /* i must be internal_relabeled */
                        internal_relabel (i);

                        if (i->d == n) break;
                        if ((l -> firstActive == sentinelNode) && 
                                        (l -> firstInactive == sentinelNode)
                           )
                                internal_gap (l);

                        if (i->d == n) break;
                }
                else {
                        /* i no longer active */
                        i->current = a;
                        /* put i on inactive list */
                        iAdd(l,i);
                        break;
                }
        } while (1);
}


// go from higher to lower buckets, push flow
void flow_solver::internal_wave() {

        node   *i;
        bucket  *l;

        for (l = buckets + aMax; l > buckets; l--) {
                for (i = l->firstActive; i != sentinelNode; i = l->firstActive) {
                        aRemove(l,i);

                        assert(i->excess > 0);
                        internal_discharge (i);

                }
        }
}

