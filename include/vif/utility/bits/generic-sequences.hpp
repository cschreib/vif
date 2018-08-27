#ifndef VIF_INCLUDING_GENERIC_BITS
#error this file is not meant to be included separately, include "vif/utilty/generic.hpp" instead
#endif

namespace vif {
    // Generate linearly increasing values.
    template<typename T = uint_t, typename ... Dims>
    vec<meta::dim_total<Dims...>::value,T> indgen(Dims&& ... ds) {
        vec<meta::dim_total<Dims...>::value,T> v(std::forward<Dims>(ds)...);
        std::iota(v.begin(), v.end(), static_cast<T>(0));
        return v;
    }
}
