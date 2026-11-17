/******************************************************************************
 * jet_refinement.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef JET_REFINEMENT_H
#define JET_REFINEMENT_H

#include "definitions.h"
#include "../refinement.h"

// Jet refinement (Gilbert, Madduri, Boman, Rajamanickam 2023), ultra only (config.jet_refinement):
// synchronous label propagation that also takes moves of bounded negative gain, filtered by the
// afterburner (gain recomputed under the higher ranked moves), alternated with a weak rebalancer
// whenever a block is overloaded. The best balanced partition seen is kept.
class jet_refinement : public refinement {
public:
        jet_refinement();
        virtual ~jet_refinement();

        virtual EdgeWeight perform_refinement(PartitionConfig & config,
                                              graph_access & G,
                                              complete_boundary & boundary);
};

#endif /* end of include guard: JET_REFINEMENT_H */
