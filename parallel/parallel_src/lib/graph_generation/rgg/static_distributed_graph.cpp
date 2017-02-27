/******************************************************************************
 * static_distributed_graph.cpp
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

#include "static_distributed_graph.h"  // This file's header.

#include <tr1/unordered_set>

#include <mpi.h>

#include "macros_assertions.h"
#include "protocol/all_to_all.h"
#include "static_local_graph.h"
#include "static_local_graph_properties.h"
#include "static_neighbour_vertex_properties.h"
#include "util.h"
#include "util_memory.h"

namespace distributed_graph {

StaticDistributedGraph::StaticDistributedGraph() :
  global_n_(kInvalidVertexId), global_m_(kInvalidEdgeId),
  local_n_(kInvalidVertexId), local_m_(kInvalidEdgeId),
  local_vertex_offset_(kInvalidVertexId), local_edge_offset_(kInvalidVertexId),
  vertex_splitters_(NULL), local_graph_(NULL), local_graph_properties_(NULL),
  neighbour_vertex_properties_(NULL)
{
}

void StaticDistributedGraph::Init(
    const MPI::Intracomm *communicator,
    long global_n,
    long global_m,
    const VertexId *vertex_splitters,
    const EdgeId *edge_splitters,
    const StaticLocalGraph *local_graph,
    StaticLocalGraphProperties *local_properties)
{
  // Duplicate communicator and run dummy global communication through it.
  // The reason for this is some funkiness for the first all-to-all
  // communication with OpenMPI 1.3.1 on IC1.
  communicator_ = communicator->Dup();

  // Initialize member variables as given parameters.
  global_n_ = global_n;
  global_m_ = global_m;
  vertex_splitters_ = vertex_splitters;
  edge_splitters_ = edge_splitters;
  local_graph_ = local_graph;
  local_graph_properties_ = local_properties;
  // Derived member variables.
  local_n_ = local_graph->n();
  local_m_ = local_graph->m();
  local_vertex_offset_ = local_graph->vertex_offset();
  local_edge_offset_ = local_graph->edge_offset();
  if (vertex_splitters_ != NULL)
    ASSERT_EQ( vertex_splitters_[communicator_.Get_rank()],
                 local_vertex_offset_ );
  // Initialize the neighbour vertex properties.
  // TODO(manuel): We only do this if vertex splitters were given and we do not use distributed graph in gap graph mode.  This is HACKY.
  if (vertex_splitters_ != NULL) 
    InitNeighbourVertexProperties();
}

StaticDistributedGraph::~StaticDistributedGraph()
{
  DeleteArrayIfAssigned(vertex_splitters_);
  DeleteArrayIfAssigned(edge_splitters_);
  DeleteObjectIfAssigned(local_graph_);
  DeleteObjectIfAssigned(local_graph_properties_);
  DeleteObjectIfAssigned(neighbour_vertex_properties_);
  communicator_.Free();
}

void StaticDistributedGraph::InitNeighbourVertexProperties()
{
  //ENTER_SECTION(&communicator_, INIT_NEIGHBOUR_VERTEX_PROPERTIES);
  //printf("Rank %3d -- InitNeighbourVertexProperties()\n",  // XXX
  //       communicator_.Get_rank());
  ProcessId process_count = communicator_.Get_size();

  // For all owned vertices v:  Send the information about v to the owners of
  // all neighbours of v.
  //
  // We transmit quadruples, encoding (v, degree, weight, internal edges).
  AllToAllHelper<long> all2all;
  all2all.Init(&communicator_);
  std::vector<std::tr1::unordered_set<VertexId> > sent_to(process_count);
  VertexId vertex_offset = local_graph_->vertex_offset();
  for (VertexId ul = 0; ul < local_n_; ++ul) {
    VertexId u = ul + vertex_offset;
    for (VertexId i = 0, iend = local_graph_->out_degree(u); i < iend; ++i) {
      VertexId v = local_graph_->neighbour(u, i);
      if (local_graph_->is_vertex_known(v))
        continue;  // Skip local vertices.
      ProcessId pid = LocateInSplitters(
          vertex_splitters_, vertex_splitters_ + process_count + 1, v);
      if (contains(sent_to[pid], u))
        continue;  // Skip if we already sent info about u there.
      sent_to[pid].insert(u);
      all2all.set_send_count(pid, all2all.send_count(pid) + 5);
    }
  }
  all2all.ExchangeSendCounts();
  typedef std::tr1::unordered_set<VertexId>::const_iterator SetIterator;
  for (ProcessId pid = 0; pid < process_count; ++pid) {
    for (SetIterator it = sent_to[pid].begin(), itend = sent_to[pid].end();
         it != itend; ++it) {
      all2all.PushForProcess(*it, pid);
      all2all.PushForProcess(local_graph_->out_degree(*it), pid);
      all2all.PushForProcess(local_graph_properties_->vertex_weight(*it), pid);
      all2all.PushForProcess(local_graph_properties_->internal_edges(*it), pid);
      // Compute interface edge weight sum.
      RggEdgeWeight w = 0;
      for (VertexCount i = 0, iend = local_graph_->out_degree(*it);
           i < iend; ++i) {
        EdgeId e = local_graph_->out_edge(*it, i);
        w += local_graph_properties_->edge_weight(e);
      }
      all2all.PushForProcess(w, pid);
    }
  }
  all2all.ExchangeData();

  // Finally, create the static neighbour vertex properties object.
  neighbour_vertex_properties_ = new StaticNeighbourVertexProperties();
  VertexCount k = all2all.receive_displacement(process_count);
  ASSERT_EQ( k % 5, 0 );
  neighbour_vertex_properties_->InitFromReceiveBuffer(
      all2all.ReceiveBuffer(), k);
  //LEAVE_SECTION(&communicator_, INIT_NEIGHBOUR_VERTEX_PROPERTIES);
}

void StaticDistributedGraph::AddUpWeights(VertexCount *n, EdgeCount *m) const
{
  *n = 0;
  *m = 0;

  // Compute local sums.
  const StaticLocalGraph *g = local_graph();
  const StaticLocalGraphProperties *p = local_graph_properties();
  VertexId vertex_offset = g->vertex_offset();
  for (VertexId u = vertex_offset, uend = vertex_offset + g->n();
       u < uend; ++u) {
    *n += p->vertex_weight(u);
    *m += p->internal_edges(u);
    for (long i = 0, iend = g->out_degree(u); i < iend; ++i) {
      *m += p->edge_weight(g->out_edge(u, i));
    }
  }

  // Compute global sums.
  long arr[2] = {*n, *m};
  communicator_.Allreduce(MPI::IN_PLACE, arr, 2, MPI::LONG, MPI::SUM);
  *n = arr[0];
  *m = arr[1];
}

}  // namespace distributed_graph

