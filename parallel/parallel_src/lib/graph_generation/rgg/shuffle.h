// File:   shuffle.h
// Author: Manuel Holtgrewe <holtgrewe@ira.uka.de>
//
// TODO(manuel): Move code to util_rand.h?

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

