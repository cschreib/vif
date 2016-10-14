#ifndef PHYPP_MATH_TRANSFORM_HPP
#define PHYPP_MATH_TRANSFORM_HPP

#include "phypp/core/vec.hpp"
#include "phypp/core/error.hpp"
#include "phypp/core/range.hpp"
#include "phypp/math/interpolate.hpp"

namespace phypp {
    // Perform the convolution of two 1D arrays, assuming that they are based on the same 'x' coordinate
    template<typename TypeX, typename TypeY1, typename TypeY2>
    auto convolve(const vec<1,TypeX>& x, const vec<1,TypeY1>& y1, const vec<1,TypeY2>& y2) ->
        vec<1,decltype(x[0]*y1[0]*y2[0])> {

        phypp_check(x.dims == y1.dims, "incompatible dimensions for X and Y1 "
            "(", x.dims, " vs. ", y1.dims, ")");
        phypp_check(x.dims == y2.dims, "incompatible dimensions for X and Y2 "
            "(", x.dims, " vs. ", y2.dims, ")");
        phypp_check(x.size() > 3, "convolve needs arrays of at least 3 elements to work "
            "(got ", x.size(), ")");

        vec<1,decltype(x[0]*y1[0]*y2[0])> r(x.size());
        for (uint_t i : range(x)) {
            auto dx = x.safe[i];
            if (i == 0) dx = x.safe[i+1] - dx;
            else dx -= x.safe[i-1];

            auto tmp = interpolate(y2, x + x.safe[i], x);
            r += y1.safe[i]*tmp*dx;
        }

        return r;
    }
}

#endif
