/******************************************************************************
 * flow_macros.h 
 *
 * based on HiPr  -- Copyright 1995, 2000 by IG Systems, Inc., igsys@eclipse.net 
 * "Maximum flow - highest lavel push-relabel algorithm"
 * Copyright for this release gpl 3.0 granted by Andrew Goldberg
 *
 * Comment: modified only the macros are left, since all the functions moved to a class named flow_solver
 *
 ******************************************************************************
 * Copyright (modifcations) (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
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


#ifndef HI_PR_5AH0M9HT
#define HI_PR_5AH0M9HT

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "values.h"
#include "types.h"          /* type definitions */
#include "timer.h"          /* timing routine */


#define GLOB_UPDT_FREQ 0.5
#define ALPHA 6
#define BETA 12

#define WHITE 0
#define GREY 1
#define BLACK 2


/* macros */

#define forAllNodes(i) for ( i = nodes; i != sentinelNode; i++ )
#define forAllArcs(i,a) for (a = i->first, stopA = (i+1)->first; a != stopA; a++)

#define nNode( i ) ( (i) - nodes + nMin )
#define nArc( a )  ( ( a == NULL )? -1 : (a) - arcs )

#define min( a, b ) ( ( (a) < (b) ) ? a : b )

#define createEdge()\
{\
                arc_first[tail + 1] ++; \
                arc_first[head + 1] ++;\
\
                /* storing information about the arc */\
                arc_tail[pos_current]        = tail;\
                arc_tail[pos_current+1]      = head;\
                arc_current       -> head    = nodes + head;\
                arc_current       -> resCap    = cap;\
                arc_current       -> rev  = arc_current + 1;\
                ( arc_current + 1 ) -> head    = nodes + tail;\
                ( arc_current + 1 ) -> resCap    = 0;\
                ( arc_current + 1 ) -> rev  = arc_current;\
\
                /* searching minimumu and maximum node */\
                if ( head < node_min ) node_min = head;\
                if ( tail < node_min ) node_min = tail;\
                if ( head > node_max ) node_max = head;\
                if ( tail > node_max ) node_max = tail;\
\
                arc_current += 2;\
                pos_current += 2;\
\
}\



#define aAdd(l,i)\
{\
  i->bNext = l->firstActive;\
  l->firstActive = i;\
  i_dist = i->d;\
  if (i_dist < aMin)\
    aMin = i_dist;\
  if (i_dist > aMax)\
    aMax = i_dist;\
  if (dMax < aMax)\
    dMax = aMax;\
}

/* i must be the first element */
#define aRemove(l,i)\
{\
  l->firstActive = i->bNext;\
}

#define iAdd(l,i)\
{\
  i_next = l->firstInactive;\
  i->bNext = i_next;\
  i->bPrev = sentinelNode;\
  i_next->bPrev = i;\
  l->firstInactive = i;\
}

#define iDelete(l,i)\
{\
  i_next = i->bNext;\
  if (l->firstInactive == i) {\
    l->firstInactive = i_next;\
    i_next->bPrev = sentinelNode;\
  }\
  else {\
    i_prev = i->bPrev;\
    i_prev->bNext = i_next;\
    i_next->bPrev = i_prev;\
  }\
}



#endif /* end of include guard: HI_PR_5AH0M9HT */
