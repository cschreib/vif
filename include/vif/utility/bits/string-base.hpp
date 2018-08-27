#ifndef PHYPP_INCLUDING_STRING_BITS
#error this file is not meant to be included separately, include "phypp/utilty/string.hpp" instead
#endif

namespace phypp {
    inline uint_t length(const std::string& s) {
        return s.size();
    }

    inline bool empty(const std::string& s) {
        return s.empty();
    }

    inline std::string keep_first(std::string s, uint_t n = 1) {
        if (s.size() > n) {
            s.erase(n);
        }

        return s;
    }

    inline std::string keep_last(std::string s, uint_t n = 1) {
        if (s.size() > n) {
            s.erase(0, s.size()-n);
        }

        return s;
    }

    inline uint_t distance(const std::string& t, const std::string& u) {
        uint_t n, d;
        if (t.size() < u.size()) {
            n = t.size();
            d = u.size() - t.size();
        } else {
            n = u.size();
            d = t.size() - u.size();
        }

        for (uint_t i : range(n)) {
            if (t[i] != u[i]) ++d;
        }

        return d;
    }

    PHYPP_VECTORIZE(length)
    PHYPP_VECTORIZE(empty)
    PHYPP_VECTORIZE(distance)
    PHYPP_VECTORIZE(keep_first)
    PHYPP_VECTORIZE(keep_last)
}
