#ifndef PHYPP_UTILITY_IDL_HPP
#define PHYPP_UTILITY_IDL_HPP

#include "phypp/core/vec.hpp"
#include "phypp/core/meta.hpp"

namespace phypp {
    // Create vectors a la IDL.
    template<typename T, typename ... Dims>
    vec<meta::dim_total<Dims...>::value, T> arr(Dims&& ... ds) {
        return vec<meta::dim_total<Dims...>::value, T>(std::forward<Dims>(ds)...);
    }

    template<typename ... Dims>
    auto fltarr(Dims ... ds) -> decltype(arr<float>(ds...)) {
        return arr<float>(ds...);
    }

    template<typename ... Dims>
    auto dblarr(Dims ... ds) -> decltype(arr<double>(ds...)) {
        return arr<double>(ds...);
    }

    template<typename ... Dims>
    auto intarr(Dims ... ds) -> decltype(arr<int_t>(ds...)) {
        return arr<int_t>(ds...);
    }

    template<typename ... Dims>
    auto uintarr(Dims ... ds) -> decltype(arr<uint_t>(ds...)) {
        return arr<uint_t>(ds...);
    }

    template<typename ... Dims>
    auto strarr(Dims ... ds) -> decltype(arr<std::string>(ds...)) {
        return arr<std::string>(ds...);
    }

    template<typename ... Dims>
    auto bytarr(Dims ... ds) -> decltype(arr<char>(ds...)) {
        return arr<char>(ds...);
    }

    template<typename ... Dims>
    auto boolarr(Dims ... ds) -> decltype(arr<bool>(ds...)) {
        return arr<bool>(ds...);
    }

    // Count the total number of elements in a vector.
    template<std::size_t Dim, typename T>
    uint_t n_elements(const vec<Dim,T>& v) {
        return v.data.size();
    }

    template<typename T, typename enable = typename std::enable_if<!meta::is_vec<T>::value>::type>
    uint_t n_elements(const T& v) {
        return 1;
    }
}


#endif
