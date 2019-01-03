/******************************************************************************
 * queue_selection_strategie.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef QUEUESELECTIONSTRATEGIE_H_
#define QUEUESELECTIONSTRATEGIE_H_

class queue_selection_strategy {
        public:
		queue_selection_strategy(PartitionConfig & config) : m_config ( config ) {};
		virtual ~queue_selection_strategy()  {};
                virtual void selectQueue(int lhs_part_weight, int rhs_part_weight, 
                                PartitionID lhs, PartitionID rhs, 
                                PartitionID & from, PartitionID & to,
                                refinement_pq * lhs_queue, refinement_pq * rhs_queue, 
                                refinement_pq** from_queue, refinement_pq** to_queue) = 0;
	protected:
		PartitionConfig m_config;

};


class queue_selection_diffusion : public queue_selection_strategy {
        public:
		queue_selection_diffusion(PartitionConfig & config) : queue_selection_strategy(config) {};
                inline void selectQueue(int lhs_part_weight, int rhs_part_weight, 
                                PartitionID lhs, PartitionID rhs, 
                                PartitionID & from, PartitionID & to,
                                refinement_pq * lhs_queue, refinement_pq * rhs_queue, 
                                refinement_pq** from_queue, refinement_pq** to_queue ) {
                        if (lhs_part_weight > rhs_part_weight) {
                                *from_queue = lhs_queue;
                                *to_queue   = rhs_queue;
                                from        = lhs;
                                to          = rhs;
                        } else {
                                *from_queue = rhs_queue;
                                *to_queue   = lhs_queue;
                                from        = rhs;
                                to          = lhs;
                        }
                }
};

class queue_selection_topgain : public queue_selection_strategy {
        public:
		queue_selection_topgain(PartitionConfig & config) : queue_selection_strategy(config) {};
                inline void selectQueue(int lhs_part_weight, int rhs_part_weight, 
                                PartitionID lhs, PartitionID rhs, 
                                PartitionID & from, PartitionID & to,
                                refinement_pq * lhs_queue, refinement_pq * rhs_queue, 
                                refinement_pq** from_queue, refinement_pq** to_queue ){

                        if( lhs_queue->empty() ) {
                                *from_queue = rhs_queue;
                                *to_queue   = lhs_queue;
                                from        = rhs;
                                to          = lhs;
                                return;
                        }
                        if( rhs_queue->empty() ) {
                                *from_queue = lhs_queue;
                                *to_queue   = rhs_queue;
                                from        = lhs;
                                to          = rhs;
                                return;
                        }

                        Gain lhsGain = lhs_queue->maxValue();
                        Gain rhsGain = rhs_queue->maxValue();

                        if(lhsGain > rhsGain){
                                *from_queue = lhs_queue;
                                *to_queue   = rhs_queue;
                                from        = lhs;
                                to          = rhs;
                        } else {
                                *from_queue = rhs_queue;
                                *to_queue   = lhs_queue;
                                from        = rhs;
                                to          = lhs;
                        }
                }
};

class queue_selection_topgain_diffusion : public queue_selection_strategy {
        public:
	  queue_selection_topgain_diffusion(PartitionConfig & config) : queue_selection_strategy(config) {  
                  qdiff = new queue_selection_diffusion(m_config);
          };

	  ~queue_selection_topgain_diffusion() {  
                  delete qdiff;
          };

          inline void selectQueue(int lhs_part_weight, int rhs_part_weight, 
                                PartitionID lhs, PartitionID rhs, 
                                PartitionID & from, PartitionID & to,
                                refinement_pq * lhs_queue, refinement_pq * rhs_queue, 
                                refinement_pq** from_queue, refinement_pq** to_queue ) {

                        if( lhs_queue->empty() ) {
                                *from_queue = rhs_queue;
                                *to_queue   = lhs_queue;
                                from        = rhs;
                                to          = lhs;
                                return;
                        }
                        if( rhs_queue->empty() ) {
                                *from_queue = lhs_queue;
                                *to_queue   = rhs_queue;
                                from        = lhs;
                                to          = rhs;
                                return;
                        }


                        Gain lhsGain = lhs_queue->maxValue();
                        Gain rhsGain = rhs_queue->maxValue();

                        if (lhsGain == rhsGain) {
                                qdiff->selectQueue(lhs_part_weight, rhs_part_weight, 
                                                   lhs, rhs, 
                                                   from, to,
                                                   lhs_queue, rhs_queue, 
                                                   from_queue, to_queue);
                                
                                return;
                        }
                        if(lhsGain > rhsGain){
                                *from_queue = lhs_queue;
                                *to_queue   = rhs_queue;
                                from        = lhs;
                                to          = rhs;
                        } else {
                                *from_queue = rhs_queue;
                                *to_queue   = lhs_queue;
                                from        = rhs;
                                to          = lhs;
                        }
                }
        private:
                queue_selection_strategy* qdiff;
};

class queue_selection_diffusion_block_targets : public queue_selection_strategy {
        public:
		queue_selection_diffusion_block_targets(PartitionConfig & config) : queue_selection_strategy(config) {
                        qdiff = new queue_selection_topgain_diffusion(config);
                };

		virtual ~queue_selection_diffusion_block_targets()  {
                        delete qdiff;
                }

                inline void selectQueue(int lhs_part_weight, int rhs_part_weight, 
                                PartitionID lhs, PartitionID rhs, 
                                PartitionID & from, PartitionID & to,
                                refinement_pq * lhs_queue, refinement_pq * rhs_queue, 
                                refinement_pq** from_queue, refinement_pq** to_queue ) {
			int lhs_overload = std::max( lhs_part_weight - m_config.target_weights[0],0);
			int rhs_overload = std::max( rhs_part_weight - m_config.target_weights[1],0);
                        if( lhs_overload == 0 && rhs_overload == 0) {
                                qdiff->selectQueue(lhs_part_weight, rhs_part_weight, 
                                                lhs, rhs, 
                                                from, to,
                                                lhs_queue, rhs_queue, 
                                                from_queue, to_queue);
                        } else {
                                if (lhs_overload > rhs_overload) {
                                        *from_queue = lhs_queue;
                                        *to_queue   = rhs_queue;
                                        from        = lhs;
                                        to          = rhs;
                                } else {
                                        *from_queue = rhs_queue;
                                        *to_queue   = lhs_queue;
                                        from        = rhs;
                                        to          = lhs;
                                }
                        }

                }

        private:
                queue_selection_strategy* qdiff;
};
#endif
