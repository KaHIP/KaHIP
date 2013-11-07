/******************************************************************************
 * flow_solver.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * based on HiPr  -- Copyright 1995, 2000 by IG Systems, Inc., igsys@eclipse.net 
 * "Maximum flow - highest lavel push-relabel algorithm"
 * Copyright for this release under gpl 3.0 granted by Andrew Goldberg
 *
 * Comment: we used the implementation of hipr and put them into a class structure here
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

#ifndef FLOW_SOLVER_4P49OMM
#define FLOW_SOLVER_4P49OMM

#include "data_structure/graph_access.h"
#include "partition_config.h"
#include "types.h"

class flow_solver {
        public:
                flow_solver( );
                virtual ~flow_solver();

                //*************************************************************************************************
                //code copied from hi_pr. make this local variables so that we can run the flow code multiple times
                //*************************************************************************************************
                void internal_stage_one ( );
                void internal_stage_two ( );

                void internal_global_update();
                void internal_check_max();
                void internal_init( );
                int  internal_allocDS( );
                void internal_wave(); 
                void internal_discharge(node* i);
                long internal_relabel(node *i);
                int  internal_gap(bucket* emptyB);

                long   n;                    /* number of nodes */
                long   m;                    /* number of arcs */
                long   nm;                   /* n + ALPHA * m */
                long   nMin;                 /* smallest node id */
                node   *nodes;               /*[> array of nodes <]*/
                node   *free_nodes;          /*[> array of nodes <]*/
                arc    *arcs;                /* array of arcs */
                bucket *buckets;             /* array of buckets */
                cType  *cap;                 /* array of capacities */
                node   *source;              /* source node pointer */
                node   *sink;                /* sink node pointer */
                long   dMax;                 /* maximum label */
                long   aMax;                 /* maximum actie node label */
                long   aMin;                 /* minimum active node label */
                double flow;                 /* flow value */
                long pushCnt;                /* number of pushes */
                long relabelCnt;             /* number of relabels */
                long updateCnt;              /* number of updates */
                long gapCnt;                 /* number of gaps */
                long gNodeCnt;               /* number of nodes after gap */  
                float t, t2;                 /* for saving times */
                node   *sentinelNode;        /* end of the node list marker */
                arc *stopA;                  /* used in forAllArcs */
                long workSinceUpdate;        /* the number of arc scans since last update */
                float globUpdtFreq;          /* global update frequency */

                long i_dist;
                node *i_next, *i_prev;
};




#endif /* end of include guard: FLOW_SOLVER_4P49OMM */

