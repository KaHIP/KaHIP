/******************************************************************************
 * triple.h
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

#ifndef TRIPLE_H
#define TRIPLE_H

namespace distributed_graph {

// A triple, similar to std::pair and easier to use than std::tr1::tuple<...>.
template <typename T1, typename T2, typename T3>
struct Triple
{
  // Default constructor.
  inline
  Triple() : first(0), second(0), third(0)
  {}

  // Construct triple from three values.
  inline
  Triple(const T1 &v1, const T2 &v2, const T3 &v3) :
    first(v1), second(v2), third(v3)
  {}

  // Lexicographic, strict weak less-than ordering.
  inline
  bool operator<(const Triple &other) const
  {
    return (first < other.first) or
      (first == other.first and
       (second < other.second or
        (second == other.second and third < other.third)));
  }

  // Lexicographic, strict weak greater-than ordering.
  inline
  bool operator>(const Triple &other) const
  {
    return (first > other.first) or
      (first == other.first and
       (second > other.second or
        (second == other.second and third > other.third)));
  }

  // First entry in the triple.
  T1 first;

  // Second entry in the triple.
  T2 second;

  // Third entry in the triple.
  T3 third;
};


// Similar to std::make_pair but for Triple.
//
// Args:
//  v1  Value for triple's first entry.
//  v2  Value for triple's second entry.
//  v3  Value for triple's third entry.
//
// Returns:
//  A a new Triple object with the given entries.
template <typename T1, typename T2, typename T3>
inline
Triple<T1, T2, T3>
MakeTriple(const T1 &v1, const T2 &v2, const T3 &v3)
{ return Triple<T1, T2, T3>(v1, v2, v3); }

}  // namespace_distributed_graph

#endif // ifndef TRIPLE_H

