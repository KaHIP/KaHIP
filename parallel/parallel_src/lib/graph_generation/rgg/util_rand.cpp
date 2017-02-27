// File:   util_rand.cpp
// Author: Manuel Holtgrewe <holtgrewe@ira.uka.de>
//
// Implementation of the PRNG-related code.

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

