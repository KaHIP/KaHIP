/******************************************************************************
 * util_rand.h
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

#ifndef UTIL_RAND_H
#define UTIL_RAND_H

#include <tr1/random>

#include <mpi.h>

namespace distributed_graph {

// Seed the given MT instance in a distributed fashion.
//
// The root process creates a PRNG and creates p random numbers.  These are
// scattered to the other processes and used to seed *mt.
//
// This function is instrumented as the section SEED_PRNG.
//
// Args:
// 
//   comm   MPI Communicator to use.
//   seed   The seed to use for generating the distributed seeds.
//   mt     The Mersenne Twister objects to use for seeding.
void SeedMersenneTwisterDistributedly(const MPI::Intracomm *comm,
                                      const int seed,
                                      std::tr1::mt19937 *mt);


// One-line function to compute a uniform integral with PRNG e.
//
// The random number will be between (including both) a and b.
//
// Args:
//
//   a  First possible value for the result.
//   b  Last possible value for the result.
//   e  PRNG engine to use for the random number generation.
template <typename Engine, typename Integral>
inline
Integral UniformInt(const Integral a, const Integral b, Engine *e)
{
  std::tr1::uniform_int<long> dist(a, b);
  return dist(*e);
}

}  // distributed_graph

#endif // ifndef UTIL_RAND_H

