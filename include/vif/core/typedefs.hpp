#ifndef PHYPP_CORE_TYPEDEF_HPP
#define PHYPP_CORE_TYPEDEF_HPP

#include <cstddef>

namespace phypp {
    // Type shortcuts
    using int_t = std::ptrdiff_t;
    using uint_t = std::size_t;

    // npos is the largest positive integer, used for signalling the end of a sequence, 
    // or referencing an invalid position in a sequence (e.g., for a non-existing element) 
    static const uint_t npos = uint_t(-1);
}

#endif
