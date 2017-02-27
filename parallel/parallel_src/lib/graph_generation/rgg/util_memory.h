/******************************************************************************
 * util_memory.h
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

