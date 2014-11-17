/******************************************************************************
 * diversifyer.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
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
