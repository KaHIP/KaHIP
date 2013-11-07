/******************************************************************************
 * timer.cpp 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * based on HiPr  -- Copyright 1995, 2000 by IG Systems, Inc., igsys@eclipse.net 
 * "Maximum flow - highest lavel push-relabel algorithm"
 * Copyright for this release gpl 3.0 granted by Andrew Goldberg
 *
 * comment: we used the implementation of hipr and put them into a header file here
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

#include "timer.h"

float timer ()
{
  struct rusage r;

  getrusage(0, &r);
  return (float)(r.ru_utime.tv_sec+r.ru_utime.tv_usec/(float)1000000);
}



