/******************************************************************************
 * types.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/


#ifndef TYPES_FIK18Y5X
#define TYPES_FIK18Y5X

#ifdef EXCESS_TYPE_LONG
typedef unsigned long excessType;
#else
typedef unsigned long long int excessType; /* change to double if not supported */
#endif

typedef long cType;
//typedef unsigned long cType;

typedef  /* arc */
   struct arcSt
{
   cType           resCap;          /* residual capasity */
   struct nodeSt   *head;           /* arc head */
   struct arcSt    *rev;            /* reverse arc */
}
  arc;

typedef  /* node */
   struct nodeSt
{
   arc             *first;           /* first outgoing arc */
   arc             *current;         /* current outgoing arc */
   excessType      excess;           /* excess at the node 
				        change to double if needed */
   long            d;                /* distance label */
   struct nodeSt   *bNext;           /* next node in bucket */
   struct nodeSt   *bPrev;           /* previous node in bucket */
} node;


typedef /* bucket */
   struct bucketSt
{
  node             *firstActive;      /* first node with positive excess */
  node             *firstInactive;    /* first node with zero excess */
} bucket;

#endif /* end of include guard: TYPES_FIK18Y5X */
