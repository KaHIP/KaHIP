/******************************************************************************
 * util.cpp
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


#include "util.h"  // This file's header.

#include <cmath>
#include <vector>


namespace distributed_graph {


void ComputeSplitters(long n, long k, long **splitters)
{
  // Number of vertices per process.
  long n_slice = n / k;
  // Vertices whose pid is < bonus_pid get a bonus element!
  long bonus_pid = n % k;
  *splitters = new long[k + 1];
  (*splitters)[0] = 0;
  for (long i = 1; i < k; ++i)
    (*splitters)[i] = (*splitters)[i - 1] + n_slice + (i < bonus_pid);
  (*splitters)[k] = n;
}


void
HsvToRgb(int h, double s, double v,
         double *r, double *g, double *b)
{
  double H, S, V, R, G, B;
  double p1, p2, p3;
  double f;
  int i;

  if (s < 0) s = 0;
  if (v < 0) v = 0;
  if (s > 1) s = 1;
  if (v > 1) v = 1;

  S = s; V = v;
  H = (h % 360) / 60.0;
  i = H;
  f = H - i;
  p1 = V * (1 - S);
  p2 = V * (1 - (S * f));
  p3 = V * (1 - (S * (1 - f)));
  if      (i == 0) { R = V;  G = p3; B = p1; }
  else if (i == 1) { R = p2; G = V;  B = p1; }
  else if (i == 2) { R = p1; G = V;  B = p3; }
  else if (i == 3) { R = p1; G = p2; B = V;  }
  else if (i == 4) { R = p3; G = p1; B = V;  }
  else             { R = V;  G = p1; B = p2; }
  *r = R;
  *g = G;
  *b = B;
}

//void RunDummyCollectiveCommunications(const MPI::Intracomm &comm)
//{
  //ENTER_SECTION(&comm, DUMMY_COLLECTIVE_COMM);
  //// Run Broadcast
  //{
    //int x = comm.Get_rank();
    //comm.Bcast(&x, 1, MPI::INTEGER, 0);
  //}
  //// Run Allgather.
  //{
    //int x = comm.Get_rank();
    //std::vector<int> rcv(comm.Get_size());
    //comm.Allgather(&x, 1, MPI::INTEGER, &rcv[0], 1, MPI::INTEGER);
  //}
  //// Run Alltoallv.
  //{
    //std::vector<int> snd(comm.Get_size());
    //std::vector<int> rcv(comm.Get_size());
    //std::vector<int> scounts(comm.Get_size(), 1);
    //std::vector<int> rcounts(comm.Get_size(), 1);
    //std::vector<int> sdispls(comm.Get_size());
    //std::vector<int> rdispls(comm.Get_size());
    //for (int i = 0, iend = sdispls.size(); i < iend; ++i) {
      //sdispls[i] = rdispls[i] = i;
    //}
    //comm.Alltoallv(&snd[0], &scounts[0], &sdispls[0], MPI::INTEGER,
                   //&rcv[0], &rcounts[0], &rdispls[0], MPI::INTEGER);
  //}
  //LEAVE_SECTION(&comm, DUMMY_COLLECTIVE_COMM);
//}

}  // namespace distributed_graph

