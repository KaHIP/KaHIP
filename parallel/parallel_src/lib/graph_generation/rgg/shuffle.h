/******************************************************************************
 * shuffle.h
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

#ifndef SHUFFLE_H
#define SHUFFLE_H

#include <tr1/random>

namespace distributed_graph {

// A simple shuffle implementation.
template <typename Iter>
void Shuffle(Iter begin, Iter end, std::tr1::mt19937 *mt)
{
  if (begin == end) return;

  for (Iter current = begin; current != end; ++current) {
    size_t n = end - current;
    std::tr1::uniform_int<size_t> dist(0, n - 1);
    std::swap(*current, *(current + dist(*mt)));
  }
}

// A simple shuffle implementation for two arrays which have to have the same length
template <typename Iter, typename IterScnd>
void Shuffle(Iter begin, Iter end, IterScnd beginScnd, IterScnd endScnd, std::tr1::mt19937 *mt)
{
  ASSERT_EQ(end - begin, endScnd - beginScnd);
  if((end - begin) != (endScnd - beginScnd)) return;
  if (begin == end) return;
   
  IterScnd currentScnd = beginScnd;

  for (Iter current = begin; current != end && currentScnd != endScnd; ++current, ++currentScnd) {
    size_t n = end - current;
    std::tr1::uniform_int<size_t> dist(0, n - 1);
    std::swap(*current, *(current + dist(*mt)));
    std::swap(*currentScnd, *(currentScnd + dist(*mt)));
  }
}


}  // namespace distributed_graph

#endif // ifndef SHUFFLE_H

