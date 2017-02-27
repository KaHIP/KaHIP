/******************************************************************************
 * helpers.h
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


#ifndef HELPERS_ZUTE7MAJ
#define HELPERS_ZUTE7MAJ

#include <algorithm>
#include <string>
#include <fstream>

class helpers {
public:
        helpers() {};
        virtual ~helpers() {};

        template<typename vectortype, class Compare, class Equal> 
        void filter_duplicates( std::vector< vectortype > & input, Compare comparator_function, Equal equal_function);
};

template<typename vectortype, class Compare, class Equal> 
void helpers::filter_duplicates( std::vector< vectortype > & input, Compare comparator_function, Equal equal_function) {
        std::sort(input.begin(), input.end(), comparator_function);

        // filter duplicates / parallel edges
        typename std::vector< vectortype >::iterator it;
        it = std::unique(input.begin(), input.end(), equal_function);

        input.resize(std::distance(input.begin(),it));

        //shrink-to-fit
        std::vector< vectortype >(input).swap(input);
}

inline bool hasEnding (std::string const &string, std::string const &ending)
{
        if (string.length() >= ending.length()) {
                return (0 == string.compare (string.length() - ending.length(), ending.length(), ending));
        } else {
                return false;
        }
}

inline bool file_exists(const std::string& name) {
        std::ifstream f(name.c_str());
        if (f.good()) {
                f.close();
                return true;
        } else {
                f.close();
                return false;
        }   
}

#endif /* end of include guard: HELPERS_ZUTE7MAJ */
