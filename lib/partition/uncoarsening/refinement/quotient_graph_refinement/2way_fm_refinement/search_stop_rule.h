/******************************************************************************
 * search_stop_rule.h 
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
