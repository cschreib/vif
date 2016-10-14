#ifndef PHYPP_MATH_RANDOM_HPP
#define PHYPP_MATH_RANDOM_HPP

#include <random>
#include "phypp/core/vec.hpp"
#include "phypp/core/error.hpp"
#include "phypp/math/base.hpp"

namespace phypp {
    using seed_t = std::mt19937;

    template<typename T>
    seed_t make_seed(T seed) {
        return std::mt19937(seed);
    }

    template<typename T>
    double randomn(T& seed) {
        std::normal_distribution<double> distribution(0.0, 1.0);
        return distribution(seed);
    }

    template<typename T, typename ... Args>
    vec<meta::dim_total<Args...>::value,double> randomn(T& seed, Args&& ... args) {
        vec<meta::dim_total<Args...>::value,double> v(std::forward<Args>(args)...);
        std::normal_distribution<double> distribution(0.0, 1.0);
        for (uint_t i : range(v)) {
            v.safe[i] = distribution(seed);
        }

        return v;
    }

    template<typename T>
    double randomu(T& seed) {
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        return distribution(seed);
    }

    template<typename T, typename ... Args>
    vec<meta::dim_total<Args...>::value,double> randomu(T& seed, Args&& ... args) {
        vec<meta::dim_total<Args...>::value,double> v(std::forward<Args>(args)...);
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        for (uint_t i : range(v)) {
            v.safe[i] = distribution(seed);
        }

        return v;
    }

    template<typename T, typename TMi, typename TMa>
    auto randomi(T& seed, TMi mi, TMa ma) -> decltype(mi + ma) {
        using rtype = decltype(mi + ma);
        std::uniform_int_distribution<rtype> distribution(mi, ma);
        return distribution(seed);
    }

    template<typename T, typename TMi, typename TMa, typename ... Args>
    auto randomi(T& seed, TMi mi, TMa ma, Args&& ... args) ->
        vec<meta::dim_total<Args...>::value,decltype(mi+ma)> {
        auto v = randomu(seed, std::forward<Args>(args)...);
        using rtype = decltype(mi + ma);
        return vec<meta::dim_total<Args...>::value,rtype>(v*(ma + 1 - mi) + mi);
    }

    template<typename T, typename TypeX, typename TypeY, typename ... Args>
    meta::rtype_t<TypeX> random_pdf(T& seed, const vec<1,TypeX>& px, const vec<1,TypeY>& py) {
        // TODO: (feature) make an alternative version for integers using std::discrete_distribution.

        using rtype = meta::rtype_t<TypeX>;
        std::piecewise_linear_distribution<rtype> distribution(px.begin(), px.end(), py.begin());
        return distribution(seed);
    }

    template<typename T, typename TypeX, typename TypeY, typename ... Args>
    vec<meta::dim_total<Args...>::value,meta::rtype_t<TypeX>> random_pdf(T& seed, const vec<1,TypeX>& px,
        const vec<1,TypeY>& py, Args&& ... args) {

        // TODO: (feature) make an alternative version for integers using std::discrete_distribution.

        using rtype = meta::rtype_t<TypeX>;
        vec<meta::dim_total<Args...>::value,rtype> v(std::forward<Args>(args)...);
        std::piecewise_linear_distribution<rtype> distribution(px.begin(), px.end(), py.begin());
        for (uint_t i : range(v)) {
            v.safe[i] = distribution(seed);
        }

        return v;
    }

    template<typename TSeed>
    bool random_coin(TSeed& seed, double prob) {
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        return distribution(seed) <= prob;
    }

    template<typename TSeed, uint_t D, typename T>
    vec<D,bool> random_coin(TSeed& seed, const vec<D,T>& prob) {
        vec<D,bool> v(prob.dims);

        std::uniform_real_distribution<double> distribution(0.0, 1.0);

        for (uint_t i : range(v)) {
            v.safe[i] = distribution(seed) <= prob.safe[i];
        }

        return v;
    }

    template<typename T, typename ... Args>
    vec<meta::dim_total<Args...>::value,bool> random_coin(T& seed, double prob, Args&& ... args) {
        vec<meta::dim_total<Args...>::value,bool> v(std::forward<Args>(args)...);
        std::uniform_real_distribution<double> distribution(0.0, 1.0);

        for (uint_t i : range(v)) {
            v.safe[i] = distribution(seed) <= prob;
        }

        return v;
    }

    template<std::size_t Dim, typename Type, typename T>
    void inplace_shuffle(T& seed, vec<Dim,Type>& v) {
        std::shuffle(v.begin(), v.end(), seed);
    }

    template<std::size_t Dim, typename Type, typename T>
    vec<Dim,Type> shuffle(T& seed, vec<Dim,Type> v) {
        inplace_shuffle(seed, v);
        return v;
    }
}

#endif
