/******************************************************************************
 * flow_macros.h 
 *
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
