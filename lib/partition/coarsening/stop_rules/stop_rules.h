/******************************************************************************
 * stop_rules.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef STOP_RULES_SZ45JQS6
#define STOP_RULES_SZ45JQS6

#include <math.h>

#include "partition_config.h"

class stop_rule {
        public:
                stop_rule() {};
                virtual ~stop_rule() {};
                virtual bool stop( NodeID number_of_finer_vertices, NodeID number_of_coarser_vertices ) = 0;
};

class separator_simple_stop_rule : public stop_rule {
        public:
                separator_simple_stop_rule(PartitionConfig & config, NodeID number_of_nodes) {
	                num_stop = config.sep_num_vert_stop;
                        if(config.disable_max_vertex_weight_constraint) {
                                config.max_vertex_weight = config.upper_bound_partition; 
                        } else {
                                config.max_vertex_weight = (NodeWeight)(1.5*config.work_load/num_stop);
                        }
                };

                virtual ~separator_simple_stop_rule() {};
                bool stop( NodeID number_of_finer_vertices, NodeID number_of_coarser_vertices );

        private:
                NodeID num_stop;
};

inline bool separator_simple_stop_rule::stop(NodeID no_of_finer_vertices, NodeID no_of_coarser_vertices ) {
        double contraction_rate = 1.0 * no_of_finer_vertices / (double)no_of_coarser_vertices;
        
        return contraction_rate >= 1.1 && no_of_coarser_vertices >= num_stop;
}

class simple_stop_rule : public stop_rule {
        public:
                simple_stop_rule(PartitionConfig & config, NodeID number_of_nodes) {
                        double x = 60;
	                num_stop = std::max(number_of_nodes/(2.0*x*config.k), 60.0*config.k);
                        if(config.disable_max_vertex_weight_constraint) {
                                config.max_vertex_weight = config.upper_bound_partition; 
                        } else {
                                config.max_vertex_weight = (NodeWeight)(1.5*config.work_load/num_stop);
                        }
                };
                virtual ~simple_stop_rule() {};
                bool stop( NodeID number_of_finer_vertices, NodeID number_of_coarser_vertices );

        private:
                NodeID num_stop;
};

inline bool simple_stop_rule::stop(NodeID no_of_finer_vertices, NodeID no_of_coarser_vertices ) {
        double contraction_rate = 1.0 * no_of_finer_vertices / (double)no_of_coarser_vertices;
        return contraction_rate >= 1.1 && no_of_coarser_vertices >= num_stop;
}
class strong_stop_rule : public stop_rule {
        public:
                strong_stop_rule(PartitionConfig & config, NodeID number_of_nodes) {
                        num_stop = config.k;
                        config.max_vertex_weight = config.upper_bound_partition;
                };
                virtual ~strong_stop_rule() {};
                bool stop( NodeID number_of_finer_vertices, NodeID number_of_coarser_vertices );

        private:
                NodeID num_stop;
};

inline bool strong_stop_rule::stop(NodeID no_of_finer_vertices, NodeID no_of_coarser_vertices ) {
        double contraction_rate = 1.0 * no_of_finer_vertices / (double)no_of_coarser_vertices;
        return contraction_rate >= 1.1 && no_of_coarser_vertices >= num_stop;
}

class multiple_k_stop_rule : public stop_rule {
        public:
                multiple_k_stop_rule (PartitionConfig & config, NodeID number_of_nodes) {
                        num_stop = config.num_vert_stop_factor*config.k;

                        if(config.disable_max_vertex_weight_constraint) {
                                config.max_vertex_weight = config.upper_bound_partition; 
                        } else {
                                if(config.initial_partitioning) {
                                        //if we perform initial partitioning we relax this constraint
                                        config.max_vertex_weight = 1.5*((double)config.work_load)/(2*config.num_vert_stop_factor); 
                                } else {
                                        config.max_vertex_weight = (NodeWeight)(1.5*config.work_load/num_stop);
                                }
                        }

                };
                virtual ~multiple_k_stop_rule () {};
                bool stop( NodeID number_of_finer_vertices, NodeID number_of_coarser_vertices );

        private:
                NodeID num_stop;
};

inline bool multiple_k_stop_rule::stop(NodeID no_of_finer_vertices, NodeID no_of_coarser_vertices ) {
        double contraction_rate = 1.0 * no_of_finer_vertices / (double)no_of_coarser_vertices;
        return contraction_rate >= 1.1 && no_of_coarser_vertices >= num_stop;
}




#endif /* end of include guard: STOP_RULES_SZ45JQS6 */
