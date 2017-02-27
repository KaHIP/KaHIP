/******************************************************************************
 * all_to_all.h
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


#ifndef ALL_TO_ALL_H
#define ALL_TO_ALL_H

#include <algorithm>
#include <numeric>
#include <utility>
#include <vector>
#include <cstdio>

#include <mpi.h>

#include "graph_types.h"
//#include "instrumentation.h"
#include "macros_assertions.h"
#include "macros_common.h"
#include "protocol/datatype_registry.h"
#include "util.h"


namespace distributed_graph {


// Type traits template for mapping data types to MPI datatypes.
template <typename T>
struct MpiDataTypeTraits
{};


// Specialization for an identifier.
//
// Note that this is based on the assumption that identifiers are 32 bit
// integers.
template <>
struct MpiDataTypeTraits<IdentifierType>
{
  static const MPI::Datatype *mpi_type()
  { return &MPI::LONG; }
};


// Specialization for a pair of 32 bit integers.
//
// This is used for sending mappings as arrays of pairs through MPI.
template <>
struct MpiDataTypeTraits<std::pair<IdentifierType, IdentifierType> >
{
  static const MPI::Datatype *mpi_type()
  { return &DatatypeRegistry::IDENTIFIER_PAIR; }
};


// A helper class for all-to-all communication with variable-length arrays.
//
// Use this where you would normally first exchange number of data elements
// to send to each processor and then the data.
//
// Usage:
//
//   AllToAllHelper<VertexId> all2all();
//   all2all.Init(&MPI::COMM_WORLD);
//   all2all.set_send_count(0, 2);
//   all2all.set_send_count(1, 3);
//   all2all.ExchangeSendCount();  // exchange and compute offsets
//   all2all.PushForProcess(0, 1);
//   all2all.PushForProcess(0, 2);
//   all2all.PushForProcess(1, 42);
//   all2all.PushForProcess(1, 2);
//   all2all.PushForProcess(1, 9);
//   all2all.ExchangeData();
//   for (ProcessId pid = 0; pid < MPI::COMM_WORLD.Get_size(); ++pid) {
//     for (int i = all2all.receive_displacement(pid);
//          i < all2all.receive_displacement(pid + 1); ++i) {
//       VertexId v = all2all.received_element(i);
//       // do something with v
//     }
//   }
//
// TODO(manuel): Error handling?
template <typename T>
class AllToAllHelper
{
public:
  // The state of the helper.  This is a simple, linear FSM.
  enum State
  {
    INVALID,
    BUILD_SEND_COUNTS,
    BUILD_DATA,
    FINISHED
  };

  // All members will be set to 0, 0.0, NULL etc.
  AllToAllHelper() :
    state_(BUILD_SEND_COUNTS), datatype_(MpiDataTypeTraits<T>::mpi_type()),
    comm_(NULL), process_count_(0),
    send_counts_(0), receive_counts_(0)
  {}

  // Initialize the helper with the given communicator.
  void Init(const MPI::Intracomm *comm)
  {
    comm_ = comm;
    process_count_ = comm_->Get_size();
    send_counts_.resize(process_count_, 0);
    receive_counts_.resize(process_count_, 0);
  }

  // Set number of elements to send to process pid to count.
  void set_send_count(ProcessId pid, int count)
  {
    ASSERT_EQ( state_, BUILD_SEND_COUNTS );
    send_counts_[pid] = count;
  }

  // Increment the number of elements to send to count by one.
  void increment_send_count(ProcessId pid)
  { set_send_count(pid, send_count(pid) + 1); }

  // Get number of elements to send to process pid.
  int send_count(ProcessId pid) const
  {
    ASSERT_EQ( state_, BUILD_SEND_COUNTS );
    return send_counts_[pid];
  }

  // Exchange the send counts.
  //
  // Precondition:  Object must be in state BUILD_SEND_COUNTS when this method
  // is called.
  //
  // Returns the number of total sent bytes.
  //
  // Postcondition:  Object is in BUILD_DATA state.
  int ExchangeSendCounts()
  {
    //REGISTER_ALL_TO_ALL(process_count_, sizeof(int) * process_count_);
    ASSERT_EQ( state_, BUILD_SEND_COUNTS );
    // XXX Uncomment the following line to make sure no additive term is measured for initializing collective communication.
    //RunDummyCollectiveCommunications(*comm_);
    // Exchange send counts.
    comm_->Alltoall(&send_counts_[0], 1, MPI::INTEGER,
                    &receive_counts_[0], 1, MPI::INTEGER);
    // Build send and receive displacements and send offsets.
    send_displacements_.resize(process_count_ + 1);
    send_displacements_[0] = 0;
    std::partial_sum(send_counts_.begin(), send_counts_.end(),
                     send_displacements_.begin() + 1);
    send_offsets_ = send_displacements_;
    send_buffer_.resize(send_displacements_.back());
    receive_displacements_.resize(process_count_ + 1);
    receive_displacements_[0] = 0;
    std::partial_sum(receive_counts_.begin(), receive_counts_.end(),
                     receive_displacements_.begin() + 1);
    receive_buffer_.resize(receive_displacements_.back());
    // Update state.
    state_ = BUILD_DATA;
    return sizeof(int) * process_count_;
  }

  // Add a new data element to send to process pid.
  void PushForProcess(const T& data, ProcessId pid)
  {
    ASSERT_EQ( state_, BUILD_DATA );
    ASSERT_BETWEEN( 0, pid, process_count_ - 1 );
    ASSERT_BETWEEN( send_displacements_[pid], send_offsets_[pid],
                      send_displacements_[pid + 1] - 1 );
    send_buffer_[send_offsets_[pid]] = data;
    send_offsets_[pid] += 1;
  }

  // Dump the send schedule for this process.
  void DumpSendSchedule(FILE *out) const
  {
    int rank = comm_->Get_rank();
    for (int i = 0; i < process_count_; ++i) {
      fprintf(out, "%4d -> %4d: %9d\n", rank, i, send_counts_[i]);
    }
  }

  // Get a pointer to the send buffer part for process pid.
  //
  // Mark the whole space as filled.
  T *SendBufferForProcess(ProcessId pid)
  {
    send_offsets_[pid] += send_counts_[pid];
    return &send_buffer_[send_displacements_[pid]];
  }

  // Get a pointer to the receive buffer part for process pid.
  T *ReceiveBufferForProcess(ProcessId pid)
  {
    return &receive_buffer_[receive_displacements_[pid]];
  }

  // Get a const pointer to the receive buffer part for process pid.
  const T *ReceiveBufferForProcess(ProcessId pid) const
  {
    return &receive_buffer_[receive_displacements_[pid]];
  }

  // Get a pointer to the overall receive buffer.
  const T *ReceiveBuffer() const
  { return &receive_buffer_[0]; }

  // Exchange the data.
  //
  // Returns number of sent bytes.
  //
  // Precondition:  Object must be in BUILD_DATA state.
  //
  // Postcondition:  Object is in FINISHED state.
  int ExchangeData()
  {
#if INSTRUMENT
    {
      int p = 0;
      int x = 0;
      for (int i = 0, iend = send_counts_.size(); i < iend; ++i) {
        x += send_counts_[i];
        p += send_counts_[i] > 0;
      }
      REGISTER_ALL_TO_ALL(p, x);
    }
#endif  // #if INSTRUMENT
#ifndef NDEBUG
    for (ProcessId pid = 0; pid < process_count_; ++pid) {
      ASSERT_EQ( send_offsets_[pid], send_displacements_[pid + 1] );
    }
#endif
    comm_->Alltoallv(&send_buffer_[0], &send_counts_[0],
                     &send_displacements_[0], *datatype_,
                     &receive_buffer_[0], &receive_counts_[0],
                     &receive_displacements_[0], *datatype_);
    state_ = FINISHED;
    return send_displacements_[process_count_] * sizeof(int);
  }

  // Reset object into BUILD_SEND_COUNTS state, reset send counts.
  void Reset()
  {
    state_ = BUILD_SEND_COUNTS;
    std::fill(send_counts_.begin(), send_counts_.end(), 0);
  }

  // TODO(manuel): Adapt to CSR-style API? received_by_len(pid)/ received_element_by(pid, i)?

  // Return the send displacement for process pid.
  int send_displacement(ProcessId pid) const
  {
    ASSERT_BETWEEN( 0, pid, process_count_ );
    return send_displacements_[pid];
  }

  // Returns number of elements sent by the process pid.
  int send_counts(ProcessId pid) const
  {
    ASSERT_BETWEEN( 0, pid, process_count_ );
    return send_counts_[pid];
  }

  // Return the i-th element sent in total.
  const T &sent_element(int i) const
  {
    ASSERT_BETWEEN( 0, i, static_cast<int>(send_buffer_.size() - 1) );
    return send_buffer_[i];
  }

  // Return the receive displacement for process pid.
  int receive_displacement(ProcessId pid) const
  {
    ASSERT_BETWEEN( 0, pid, process_count_ );
    return receive_displacements_[pid];
  }

  // Returns number of elements received by the process pid.
  int receive_counts(ProcessId pid) const
  {
    ASSERT_BETWEEN( 0, pid, process_count_ );
    return receive_counts_[pid];
  }

  // Return the i-th element received in total.
  const T &received_element(int i) const
  {
    ASSERT_EQ( state_, FINISHED );
    ASSERT_BETWEEN( 0, i, static_cast<int>(receive_buffer_.size() - 1) );
    return receive_buffer_[i];
  }

  int received_count() const
  {
    ASSERT_EQ( state_, FINISHED );
    return receive_displacements_[process_count_];
  }

private:
  // The object's current state.
  State state_;

  // We will send arrays of this MPI datatype.
  const MPI::Datatype *datatype_;

  // The communicator we work in.
  const MPI::Intracomm *comm_;

  // The number of processes, for brevity in the functions above.
  ProcessId process_count_;

  // The send buffer.
  std::vector<T> send_buffer_;

  // The number of elements to send to each process.
  std::vector<int> send_counts_;

  // The displacements in the final send buffer, fixed after they are once
  // computed.
  std::vector<int> send_displacements_;

  // Offsets in the send buffer, adjusted in PushForProcess().
  std::vector<int> send_offsets_;

  // The receive buffer.
  std::vector<T> receive_buffer_;

  // The number of elements received by each process.
  std::vector<int> receive_counts_;

  // The displacements in the receive buffer.
  std::vector<int> receive_displacements_;

  DISALLOW_COPY_AND_ASSIGN(AllToAllHelper);
};


}  // namespace distributed_graph


#endif /* ifndef ALL_TO_ALL_H */

