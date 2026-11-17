/******************************************************************************
 * aae_timers.h
 * *
 * Per-phase wall-clock timers per process (profiling instrumentation).
 * Each process accumulates the time it spends in each phase; at the end the
 * maximum and the mean over all processes are printed by rank 0 as
 * "aae_time <phase> <max> <mean>". The timers only read the clock, so they
 * do not change any result.
 *****************************************************************************/

#ifndef AAE_TIMERS_PARHIP_H
#define AAE_TIMERS_PARHIP_H

#include <chrono>
#include <iostream>
#include <vector>
#include <mpi.h>

namespace parhip {

enum aae_phase {
        AAE_IO = 0,       // reading the graph
        AAE_COARSEN_LP,   // label propagation clustering during coarsening (incl. balance setup)
        AAE_CONTRACT,     // contraction (and block down propagation in later V-cycles)
        AAE_IP,           // initial partitioning: collect, distribute, kaffpaE_comm, assignment
        AAE_IP_KAFFPAE,   // the kaffpaE_comm call alone (part of AAE_IP)
        AAE_PROJECT,      // projection to the finer level
        AAE_REFINE_LP,    // label propagation refinement
        AAE_COMM,         // ghost node exchange inside label propagation (part of the two LP phases)
        AAE_WRITE,        // writing the partition file
        AAE_NUM_PHASES
};

inline double* aae_phase_times() {
        static double t[AAE_NUM_PHASES] = {0};
        return t;
}

inline double aae_now() {
        return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

struct aae_scope {
        int phase;
        double start;
        explicit aae_scope(int p) : phase(p), start(aae_now()) {}
        ~aae_scope() { aae_phase_times()[phase] += aae_now() - start; }
};

inline void aae_report(MPI_Comm comm) {
        static const char* names[AAE_NUM_PHASES] = {"io", "coarsen_lp", "contract", "ip", "ip_kaffpaE",
                                                    "project", "refine_lp", "comm", "write"};
        int rank, size;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &size);
        std::vector<double> mx(AAE_NUM_PHASES), sum(AAE_NUM_PHASES);
        MPI_Reduce(aae_phase_times(), &mx[0], AAE_NUM_PHASES, MPI_DOUBLE, MPI_MAX, 0, comm);
        MPI_Reduce(aae_phase_times(), &sum[0], AAE_NUM_PHASES, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (rank == 0) {
                for (int i = 0; i < AAE_NUM_PHASES; i++) {
                        std::cout << "aae_time " << names[i] << " " << mx[i] << " " << sum[i] / size << std::endl;
                }
        }
}

} // namespace parhip

using parhip::aae_scope;
using parhip::aae_report;
using parhip::AAE_IO;
using parhip::AAE_COARSEN_LP;
using parhip::AAE_CONTRACT;
using parhip::AAE_IP;
using parhip::AAE_IP_KAFFPAE;
using parhip::AAE_PROJECT;
using parhip::AAE_REFINE_LP;
using parhip::AAE_COMM;
using parhip::AAE_WRITE;

#endif
