/******************************************************************************
 * util_rand.cpp
 *
 * Source of KaHIP -- Karlsruhe High Quality Graph Partitioning 
 ******************************************************************************
 * Copyright (C) 2017 Christian Schulz 
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

#include "util_rand.h"  // This file's header.

//#include "instrumentation.h"

namespace distributed_graph {

void SeedMersenneTwisterDistributedly(const MPI::Intracomm *comm,
                                      const int seed,
                                      std::tr1::mt19937 *mt)
{
  //ENTER_SECTION(comm, SEED_PRNG);
  int process_count = comm->Get_size();
  int rank = comm->Get_rank();

  // The root process creates a PRNG and creates p random numbers.  These are
  // scattered to the other processes and used to seed *mt.
  int *global_seeds = NULL;
  if (rank == 0) {
    global_seeds = new int[process_count];
    std::tr1::mt19937 initial_mt(seed);
    for (int i = 0; i < process_count; ++i)
      global_seeds[i] = initial_mt();
  }

  int local_seed;
  comm->Scatter(global_seeds, 1, MPI::INTEGER,
                &local_seed, 1, MPI::INTEGER, 0);
  mt->seed(local_seed);

  // Cleanup.
  if (global_seeds != NULL)
    delete[] global_seeds;
  //LEAVE_SECTION(comm, SEED_PRNG);
}

}  // namespace distributed_graph

