/******************************************************************************
 * partition_accept_rule.h 
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

#ifndef PARTITION_ACCEPT_RULE_4RXUS4P9
#define PARTITION_ACCEPT_RULE_4RXUS4P9

#include "partition_config.h"
#include "random_functions.h"

class partition_accept_rule {
        public:
                partition_accept_rule( ) {};
                virtual ~partition_accept_rule() {};

                virtual bool accept_partition(PartitionConfig & config, 
                                              const EdgeWeight edge_cut, 
                                              const NodeWeight lhs_part_weight, 
                                              const NodeWeight rhs_part_weight,
                const PartitionID lhs, 
                const PartitionID rhs, 
                                              bool & rebalance ) = 0;


};

class normal_partition_accept_rule : public partition_accept_rule {
        public:
                normal_partition_accept_rule(PartitionConfig & config, 
                                             const EdgeWeight initial_cut, 
                                             const NodeWeight initial_lhs_part_weight,
                                             const NodeWeight initial_rhs_part_weight);
                virtual ~normal_partition_accept_rule() {};

                bool accept_partition(PartitionConfig & config, 
                                      const EdgeWeight edge_cut, 
                                      const NodeWeight lhs_part_weight,
                                      const NodeWeight rhs_part_weight, 
                const PartitionID lhs, 
                const PartitionID rhs, 

                                      bool & rebalance);
        private:
                EdgeWeight best_cut;
                NodeWeight cur_lhs_part_weight;
                NodeWeight cur_rhs_part_weight;
                NodeWeight difference;
};

normal_partition_accept_rule::normal_partition_accept_rule(PartitionConfig & config, 
                const EdgeWeight initial_cut, 
                const NodeWeight initial_lhs_part_weight, 
                const NodeWeight initial_rhs_part_weight) {

        best_cut            = initial_cut;
        cur_lhs_part_weight = initial_lhs_part_weight;
        cur_rhs_part_weight = initial_rhs_part_weight;
        difference          = abs((int)cur_lhs_part_weight - (int)cur_rhs_part_weight);
}

bool normal_partition_accept_rule::accept_partition(PartitionConfig & config, 
                const EdgeWeight edge_cut, 
                const NodeWeight lhs_part_weight, 
                const NodeWeight rhs_part_weight,
                const PartitionID lhs, 
                const PartitionID rhs, 
                bool & rebalance) {
        NodeWeight cur_diff            = abs((int)lhs_part_weight - (int)rhs_part_weight);
        bool better_cut_within_balance = edge_cut < best_cut;

        if(config.softrebalance) {
              better_cut_within_balance = edge_cut <= best_cut;  
        }

        better_cut_within_balance = better_cut_within_balance && 
                                    lhs_part_weight < config.upper_bound_partition 
                                 && rhs_part_weight < config.upper_bound_partition; 

        if( (better_cut_within_balance 
         || (cur_diff < difference && edge_cut == best_cut)) 
         &&  lhs_part_weight > 0 && rhs_part_weight > 0 ) {

                best_cut   = edge_cut;
                difference = cur_diff;
                rebalance  = false;
                return true; 

        } else if(rebalance) {
                if(cur_diff < difference 
                || (cur_diff <= difference && edge_cut < best_cut)) {
                        best_cut   = edge_cut;
                        difference = cur_diff;
                        return true;
                }        
        }
        return false;        
}

class ip_partition_accept_rule : public partition_accept_rule {
        public:
                ip_partition_accept_rule(PartitionConfig & config, 
                                             const EdgeWeight initial_cut, 
                                             const NodeWeight initial_lhs_part_weight,
                                             const NodeWeight initial_rhs_part_weight, 
                                             const PartitionID lhs, 
                                            const PartitionID rhs);
                virtual ~ip_partition_accept_rule() {};

                bool accept_partition(PartitionConfig & config, 
                                      const EdgeWeight edge_cut, 
                                      const NodeWeight lhs_part_weight,
                                      const NodeWeight rhs_part_weight, 
                const PartitionID lhs, 
                const PartitionID rhs, 
                                      bool & rebalance);
        private:
                EdgeWeight best_cut;
                int cur_lhs_overload;
                int cur_rhs_overload;
};

ip_partition_accept_rule::ip_partition_accept_rule(PartitionConfig & config, 
                const EdgeWeight initial_cut, 
                const NodeWeight initial_lhs_part_weight, 
                const NodeWeight initial_rhs_part_weight,
                const PartitionID lhs, 
                const PartitionID rhs) {

        best_cut            = initial_cut;
        cur_lhs_overload = std::max( (int)initial_lhs_part_weight - config.target_weights[lhs],0);
        cur_rhs_overload = std::max( (int)initial_rhs_part_weight - config.target_weights[rhs],0);
}

bool ip_partition_accept_rule::accept_partition(PartitionConfig & config, 
                const EdgeWeight edge_cut, 
                const NodeWeight lhs_part_weight, 
                const NodeWeight rhs_part_weight,
                const PartitionID lhs, 
                const PartitionID rhs, 
                bool & rebalance) {
        bool better_cut_within_balance = edge_cut <= best_cut;

        int act_lhs_overload = std::max( (int)lhs_part_weight - config.target_weights[lhs],0);
        int act_rhs_overload = std::max( (int)rhs_part_weight - config.target_weights[rhs],0);

        better_cut_within_balance = better_cut_within_balance && 
                                    act_lhs_overload == 0 && act_rhs_overload == 0; 

        if( act_rhs_overload == 0 && act_lhs_overload == 0 ) config.rebalance = false;

        if( config.rebalance ) {
                if(act_rhs_overload + act_lhs_overload < cur_lhs_overload + cur_rhs_overload
                || (act_rhs_overload + act_lhs_overload <= cur_lhs_overload + cur_rhs_overload
                    && edge_cut < best_cut)) {
                        best_cut   = edge_cut;
                        cur_lhs_overload = act_lhs_overload;
                        cur_rhs_overload = act_rhs_overload;
                        return true;
                }        
        } else {
                if( (better_cut_within_balance 
                    || (act_rhs_overload + act_lhs_overload < cur_lhs_overload + cur_rhs_overload && edge_cut == best_cut)) 
                       &&  lhs_part_weight > 0 && rhs_part_weight > 0 ) {

                        best_cut         = edge_cut;
                        cur_lhs_overload = act_lhs_overload;
                        cur_rhs_overload = act_rhs_overload;
                        return true; 

                } 
        }
        return false;        
}

#endif /* end of include guard: PARTITION_ACCEPT_RULE_4RXUS4P9 */
