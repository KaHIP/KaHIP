/******************************************************************************
 * diversifyer.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef DIVERSIFYER_AZQIF42R
#define DIVERSIFYER_AZQIF42R

#include "random_functions.h"

class diversifyer {
public:
        diversifyer() {} ;
        virtual ~diversifyer() {};

        void diversify(PartitionConfig & config) {
                //diversify edge rating:
                config.edge_rating                   = (EdgeRating)random_functions::nextInt(0, (unsigned)EXPANSIONSTAR2ALGDIST);
                config.permutation_quality           = PERMUTATION_QUALITY_GOOD;
                config.permutation_during_refinement = PERMUTATION_QUALITY_GOOD;
        }

        void diversify_kaba(PartitionConfig & config) {
                config.kaba_flip_packings             = false;
                config.kaba_packing_iterations        = random_functions::nextInt(1, 20);
                config.kaba_internal_no_aug_steps_aug = random_functions::nextInt(1, 30);
                config.kaba_unsucc_iterations         = random_functions::nextInt(1, 10);
        }
};


#endif /* end of include guard: DIVERSIFYER_AZQIF42R */
