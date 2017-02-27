// File:   triple.h
// Author: Manuel Holtgrewe <holtgrewe@ira.uka.de>
//
// A std::pair<>-like Triple class.

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

