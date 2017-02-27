/******************************************************************************
 * util.h
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

#ifndef UTIL_H
#define UTIL_H

#include <algorithm>          // lower_bound()
#include <cmath>              // ceil()
#include <cstdlib>
#include <functional>         // unary_function<,>, binary_function<,,>
#include <tr1/unordered_map>  // hash<>
#include <utility>            // pair<,>
#include <mpi.h>

#include "macros_assertions.h"

namespace distributed_graph {

inline
void AbortProgram(const char * message)
{
  printf("ABORT:  %s\n", message);
  abort();
}

// Absolute values for integers.
inline
long iabs(long x)
{ return x >= 0 ? x : -x; }


// Template that helps in building an increment operator for enums.
//
// Example:
//
// enum MyEnum { BEGIN, VALUE1 = BEGIN, VALUE2, END };
//
// inline
// MyEnum &operator++(MyEnum &v)
// { return enum_increment(v, BEGIN, END); }
//
// for (MyEnum e = BEGIN; e != END; ++e) {
//   // ...
// }
template <class Enum>
Enum & enum_increment(Enum & value, Enum begin, Enum end)
{ return value = (value == end) ? begin : Enum(value + 1); }


// Direction for comparisons.
enum CompareDirection
{
  COMPARE_LESS,
  COMPARE_GREATER
};


// Fill a sequence with incrementing values.
//
// The sequence will be filled with the values (first, first + 1, first + 2,
// first + 3, ...).
//
// Args:
//   begin  Iterator to first value of the sequence.
//   end    Iterator behind the last value of the sequence.
//   first  The first value to assign.
template <typename Iter, typename Integral>
inline
void FillWithIncrementingValues(Iter begin, Iter end, Integral first)
{
  for (Iter it = begin; it != end; ++it, ++first) {
    *it = first;
  }
}


// Functor to compare integer values by looking them up in a table/vector.
//
// Usage:
//
//   int *values = {3, 0, 1};
//   int *keys = {0, 1, 2};
//   std::sort(keys, keys + 3, CompareAsKey<int*, int>(values);
template <typename Iter, typename Integral, CompareDirection lt = COMPARE_LESS>
class CompareAsKey : public std::binary_function<long, long, bool>
{
public:
  // Initialize the comparison functor with the key sequence beginning at
  // values.
  // TODO(manuel): No default values!
  CompareAsKey(Iter values, Integral offset = 0) : values_(values), offset_(offset)
  {}

  // Compare the values in the key sequence at position left with the value
  // at position right.
  //
  // The comparison direction (< vs. >) depends on the template parameter lt.
  bool operator()(const Integral &left, const Integral &right) const
  {
    if (lt == COMPARE_LESS)
      return values_[left - offset_] < values_[right - offset_];
    else
      return values_[left - offset_] > values_[right - offset_];
  }

private:

  // Iterator to the first value in the value sequence.
  Iter values_;

  // The offset of the keys.
  Integral offset_;
};


// Compute the splitters for n elements when separating into k parts.
//
// This function creates the splitters array in *splitters.  It is your
// responsibility to free them with delete[].
void ComputeSplitters(long n, long k, long **splitters);


// "Locate" the given value in the splitters sequence.
//
// This means it returns the index of the block it would fall into.
//
// Consider the following examples.
//
//   int splitters[] = {0, 2, 5, 9};
//   LocateInSplitters(splitters, splitters + 4, 0);  // => 0
//   LocateInSplitters(splitters, splitters + 4, 1);  // => 0
//   LocateInSplitters(splitters, splitters + 4, 2);  // => 1
//   LocateInSplitters(splitters, splitters + 4, 3);  // => 1
//   LocateInSplitters(splitters, splitters + 4, 5);  // => 2
template <typename Iter, typename T>
inline
int LocateInSplitters(Iter begin, Iter end, T value)
{
  ASSERT_TRUE( begin != end );
  ASSERT_BETWEEN( *begin, value, *(end - 1) - 1 );
  Iter found = std::upper_bound(begin, end, value);
  return found - begin - 1;
}


// Convert the given HSV value into RGB.
//
// Taken from xscreensaver.
//
// TODO(manuel): Add copyright notice?
void HsvToRgb(int h, double s, double v,
              double *r, double *g, double *b);


// Hashing for pairs.
//
// Usage:
//
//   std::tr1::unordered_set<
//       std::pair<int, int>,
//       PairHash<std::pair<int, int> > myset;
template <typename T1, typename T2>
struct PairHash : public std::unary_function<std::pair<T1, T2>, size_t>
{
  size_t operator()(const ::std::pair<T1, T2> &p) const
  { return std::tr1::hash<T1>()(p.first) + std::tr1::hash<T2>()(p.second); }
};


// A "contains" operator for STL containers.
//
// Equivalent to (c.find(x) != c.end()).
template <typename Container, typename T>
inline
bool contains(const Container &c, const T &x)
{ return c.find(x) != c.end(); }

// Run some dummy collective communications.
void RunDummyCollectiveCommunications(const MPI::Intracomm &comm);

}  // namespace distributed_graph

#endif // ifndef UTIL_H

