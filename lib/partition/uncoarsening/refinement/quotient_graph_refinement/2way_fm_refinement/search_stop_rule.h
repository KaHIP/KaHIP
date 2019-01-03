/******************************************************************************
 * search_stop_rule.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef SEARCH_STOP_RULE_R20GH6IN
#define SEARCH_STOP_RULE_R20GH6IN

class stop_rule {
        public:
                stop_rule( ) {};
                virtual ~stop_rule() {};

                virtual bool search_should_stop(unsigned int min_cut_idx, 
                                                unsigned int cur_idx, 
                                                unsigned int search_limit) = 0;
};

class easy_stop_rule : public stop_rule {
        public:
                easy_stop_rule( ) {};
                virtual ~easy_stop_rule() {};

                bool search_should_stop(unsigned int min_cut_idx, 
                                        unsigned int cur_idx, 
                                        unsigned int search_limit);
};

inline bool easy_stop_rule::search_should_stop(unsigned min_cut_idx, 
                                               unsigned int cur_idx, 
                                               unsigned int search_limit) {
        return cur_idx - min_cut_idx > search_limit;
}

#endif /* end of include guard: SEARCH_STOP_RULE_R20GH6IN */
