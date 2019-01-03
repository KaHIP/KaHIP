/******************************************************************************
 * helpers.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
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
