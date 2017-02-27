/******************************************************************************
 * types.h 
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
