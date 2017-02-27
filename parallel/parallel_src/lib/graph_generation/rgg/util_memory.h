// File:   util_memory.h
// Author: Manuel Holtgrewe <holtgrewe@ira.uka.de>
//
// Memory management related utils:
//
// Functions:
//
//   DeleteObjectIfAssigned()
//   DeleteArrayIfAssigned()

#ifndef UTIL_MEMORY_H
#define UTIL_MEMORY_H

namespace distributed_graph {

// Call delete on ptr if ptr is not NULL.
//
// Arguments:
//
//   ptr  Pointer to the object to delete.
template <typename T>
inline
void DeleteObjectIfAssigned(T *ptr)
{ if (ptr != NULL) delete ptr; }


// Call delete[] on ptr if ptr is not NULL.
//
// Arguments:
//
//   ptr  Pointer to the array to delete.
template <typename T>
inline
void DeleteArrayIfAssigned(T *ptr)
{ if (ptr != NULL) delete[] ptr; }


// Call delete on all pointers in the given sequence if the value the pointer
// points to is non-NULL.
//
// Arguments:
//
//   begin  Iterator to first element in the sequence.
//   end    Iterator beyond the last element in the sequence.
template <typename Iter>
inline
void DeleteObjectsInSequenceIfAssigned(Iter begin, Iter end)
{
  for (Iter it = begin; it != end; ++it) {
    if (*it != NULL) delete *it;
  }
}

}  // namespace distributed_graph

#endif // ifndef UTIL_MEMORY_H

