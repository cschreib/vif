#ifndef PHYPP_MATH_MATH_HPP
#define PHYPP_MATH_MATH_HPP

#include <cmath>
#ifndef NO_GSL
#include <gsl/gsl_sf_bessel.h>
#endif
#include "phypp/core/vec.hpp"
#include "phypp/core/range.hpp"
#include "phypp/core/error.hpp"

namespace phypp {
    static constexpr const double dnan = std::numeric_limits<double>::quiet_NaN();
    static constexpr const float  fnan = std::numeric_limits<float>::quiet_NaN();
    static constexpr const double dinf = std::numeric_limits<double>::infinity();
    static constexpr const float  finf = std::numeric_limits<float>::infinity();
    static constexpr const double dpi = 3.14159265359;
    static constexpr const float  fpi = 3.14159265359;

    // Import some standard functions into the global namespace for convenience
    using std::sqrt;
    using std::cos;
    using std::sin;
    using std::tan;
    using std::acos;
    using std::asin;
    using std::atan;
    using std::atan2;
    using std::cosh;
    using std::sinh;
    using std::tanh;
    using std::acosh;
    using std::asinh;
    using std::atanh;
    using std::pow;
    using std::fmod;
    using std::exp;
    using std::log;
    using std::log2;
    using std::log10;
    using std::erf;
    using std::erfc;
    using std::tgamma;
    using std::ceil;
    using std::floor;
    using std::round;
    using std::abs;

    template<typename T, typename enable = typename std::enable_if<!meta::is_vec<T>::value>::type>
    bool is_finite(const T& t) {
        return std::isfinite(t);
    }

    template<std::size_t Dim, typename Type>
    vec<Dim,bool> is_finite(const vec<Dim,Type>& v) {
        vec<Dim,bool> r(v.dims);
        for (uint_t i : range(v)) {
            r.safe[i] = std::isfinite(v.safe[i]);
        }

        return r;
    }

    template<typename T, typename enable = typename std::enable_if<!meta::is_vec<T>::value>::type>
    bool is_nan(const T& t) {
        return std::isnan(t);
    }

    template<std::size_t Dim, typename Type>
    vec<Dim,bool> is_nan(const vec<Dim,Type>& v) {
        vec<Dim,bool> r(v.dims);
        for (uint_t i : range(v)) {
            r.safe[i] = std::isnan(v.safe[i]);
        }

        return r;
    }

    template<typename T>
    auto e10(const T& t) -> decltype(pow(10.0, t)) {
        return pow(10.0, t);
    }

    template<typename T, typename U, typename V,
        typename enable = typename std::enable_if<!meta::is_vec<T>::value>::type>
    T clamp(const T& t, const U& mi, const V& ma) {
        return (t < mi ? mi : (t > ma ? ma : t));
    }

    template<typename T>
    auto sign(const T& t) -> decltype(2*(t >= 0) - 1) {
        return 2*(t >= 0) - 1;
    }

    template<std::size_t Dim, typename Type, typename U,
        typename enable = typename std::enable_if<!meta::is_vec<U>::value>::type>
    auto pow(const U& u, const vec<Dim,Type>& v) -> vec<Dim, decltype(pow(u,v(0)))> {
        vec<Dim, decltype(pow(u,v(0)))> r = v;
        for (auto& t : r) {
            t = pow(u, t);
        }
        return r;
    }

    template<std::size_t Dim, typename Type, typename U>
    auto pow(const U& u, vec<Dim,Type>&& v) -> typename std::enable_if<!meta::is_vec<U>::value &&
        !std::is_pointer<Type>::value && std::is_same<decltype(pow(u,v(0))), Type>::value,
        vec<Dim,Type>>::type {
        for (auto& t : v) {
            t = pow(u, t);
        }
        return std::move(v);
    }

    template<typename T, typename enable = typename std::enable_if<!meta::is_vec<T>::value>::type>
    auto sqr(T t) -> decltype(t*t) {
        return t*t;
    }

    template<typename T, typename enable = typename std::enable_if<!meta::is_vec<T>::value>::type>
    auto invsqr(T t) -> decltype(1.0/(t*t)) {
        return 1.0/(t*t);
    }

    #define VECTORIZE(name) \
        template<std::size_t Dim, typename Type, typename ... Args> \
        auto name(const vec<Dim,Type>& v, const Args& ... args) -> \
            vec<Dim,decltype(name(v[0], args...))> { \
            using ntype = decltype(name(v[0], args...)); \
            vec<Dim,ntype> r; r.dims = v.dims; r.data.reserve(v.size()); \
            for (auto& t : v.data) { \
                r.data.push_back(name(impl::dref<Type>(t), args...)); \
            } \
            return r; \
        } \
        template<std::size_t Dim, typename Type, typename ... Args> \
        auto name(vec<Dim,Type>&& v, const Args& ... args) -> typename std::enable_if< \
            !std::is_pointer<Type>::value && std::is_same<decltype(name(v[0], args...)), Type>::value, \
            vec<Dim,Type>>::type { \
            for (auto& t : v) { \
                t = name(t, args...); \
            } \
            return std::move(v); \
        }

    #define VECTORIZE2(name) \
        template<std::size_t D, typename T1, typename T2, typename ... Args> \
        auto name(const vec<D,T1>& v1, const vec<D,T2>& v2, const Args& ... args) -> \
            vec<D,decltype(name(v1[0], v2[0], args...))> { \
            phypp_check(v1.dims == v2.dims, "incompatible dimensions between V1 and V2 (", \
                v1.dims, " vs. ", v2.dims, ")"); \
            using ntype = decltype(name(v1[0], v2[0], args...)); \
            vec<D,ntype> r; r.dims = v1.dims; r.data.reserve(v1.size()); \
            for (uint_t i : range(v1)) { \
                r.data.push_back(name(v1.safe[i], v2.safe[i], args...)); \
            } \
            return r; \
        } \
        template<std::size_t D, typename T1, typename T2, typename ... Args> \
        auto name(vec<D,T1>&& v1, const vec<D,T2>& v2, const Args& ... args) -> typename std::enable_if< \
            !std::is_pointer<T1>::value && std::is_same<decltype(name(v1[0], v2[0], args...)), T1>::value, \
            vec<D,T1>>::type { \
            phypp_check(v1.dims == v2.dims, "incompatible dimensions between V1 and V2 (", \
                v1.dims, " vs. ", v2.dims, ")"); \
            for (uint_t i : range(v1)) { \
                v1.safe[i] = name(v1.safe[i], v2.safe[i], args...); \
            } \
            return v1; \
        } \
        template<std::size_t D, typename T1, typename T2, typename ... Args> \
        auto name(const vec<D,T1>& v1, vec<D,T2>&& v2, const Args& ... args) -> typename std::enable_if< \
            !std::is_pointer<T2>::value && std::is_same<decltype(name(v1[0], v2[0], args...)), T2>::value, \
            vec<D,T2>>::type { \
            phypp_check(v1.dims == v2.dims, "incompatible dimensions between V1 and V2 (", \
                v1.dims, " vs. ", v2.dims, ")"); \
            for (uint_t i : range(v1)) { \
                v2.safe[i] = name(v1.safe[i], v2.safe[i], args...); \
            } \
            return v2; \
        } \
        template<std::size_t D, typename T1, typename T2, typename ... Args> \
        auto name(vec<D,T1>&& v1, vec<D,T2>&& v2, const Args& ... args) -> typename std::enable_if< \
            !std::is_pointer<T1>::value && std::is_same<decltype(name(v1[0], v2[0], args...)), T1>::value, \
            vec<D,T1>>::type { \
            phypp_check(v1.dims == v2.dims, "incompatible dimensions between V1 and V2 (", \
                v1.dims, " vs. ", v2.dims, ")"); \
            for (uint_t i : range(v1)) { \
                v1.safe[i] = name(v1.safe[i], v2.safe[i], args...); \
            } \
            return v1; \
        } \
        template<std::size_t D, typename T1, typename T2, typename ... Args> \
        auto name(T1 v1, const vec<D,T2>& v2, const Args& ... args) -> typename std::enable_if<!meta::is_vec<T1>::value, \
            vec<D,decltype(name(v1, v2[0], args...))>>::type { \
            using ntype = decltype(name(v1, v2[0], args...)); \
            vec<D,ntype> r; r.dims = v2.dims; r.data.reserve(v2.size()); \
            for (uint_t i : range(v2)) { \
                r.data.push_back(name(v1, v2.safe[i], args...)); \
            } \
            return r; \
        } \
        template<std::size_t D, typename T1, typename T2, typename ... Args> \
        auto name(const vec<D,T1>& v1, T2 v2, const Args& ... args) -> typename std::enable_if<!meta::is_vec<T2>::value, \
            vec<D,decltype(name(v1[0], v2, args...))>>::type { \
            using ntype = decltype(name(v1[0], v2, args...)); \
            vec<D,ntype> r; r.dims = v1.dims; r.data.reserve(v1.size()); \
            for (uint_t i : range(v1)) { \
                r.data.push_back(name(v1.safe[i], v2, args...)); \
            } \
            return r; \
        } \
        template<std::size_t D, typename T1, typename T2, typename ... Args> \
        auto name(T1 v1, vec<D,T2>&& v2, const Args& ... args) -> typename std::enable_if<!meta::is_vec<T1>::value && \
            !std::is_pointer<T2>::value && std::is_same<decltype(name(v1, v2[0], args...)), T2>::value, \
            vec<D,decltype(name(v1, v2[0], args...))>>::type { \
            for (uint_t i : range(v2)) { \
                v2.safe[i] = name(v1, v2.safe[i], args...); \
            } \
            return v2; \
        } \
        template<std::size_t D, typename T1, typename T2, typename ... Args> \
        auto name(vec<D,T1>&& v1, T2 v2, const Args& ... args) -> typename std::enable_if<!meta::is_vec<T2>::value && \
            !std::is_pointer<T1>::value && std::is_same<decltype(name(v1[0], v2, args...)), T1>::value, \
            vec<D,decltype(name(v1[0], v2, args...))>>::type { \
            for (uint_t i : range(v1)) { \
                v1.safe[i] = name(v1.safe[i], v2, args...); \
            } \
            return v1; \
        } \

    #define VECTORIZE_REN(name, orig) \
        template<std::size_t Dim, typename Type, typename ... Args> \
        auto name(const vec<Dim,Type>& v, const Args& ... args) -> \
            vec<Dim,decltype(orig(v[0], args...))> { \
            using ntype = decltype(orig(v[0], args...)); \
            vec<Dim,ntype> r; r.dims = v.dims; r.data.reserve(v.size()); \
            for (auto& t : v.data) { \
                r.data.push_back(orig(impl::dref<Type>(t), args...)); \
            } \
            return r; \
        } \
        template<std::size_t Dim, typename Type, typename ... Args> \
        auto name(vec<Dim,Type>&& v, const Args& ... args) -> typename std::enable_if< \
            !std::is_pointer<Type>::value && std::is_same<decltype(orig(v[0], args...)), Type>::value, \
            vec<Dim,Type>>::type { \
            for (auto& t : v) { \
                t = orig(t, args...); \
            } \
            return std::move(v); \
        } \
        template<typename ... Args> \
        auto name(Args&& ... args) -> decltype(orig(std::forward<Args>(args)...)) { \
            return orig(std::forward<Args>(args)...); \
        }

    VECTORIZE(sqrt);
    VECTORIZE(sqr);
    VECTORIZE(invsqr);
    VECTORIZE(pow);
    VECTORIZE(fmod);
    VECTORIZE(cos);
    VECTORIZE(sin);
    VECTORIZE(tan);
    VECTORIZE(acos);
    VECTORIZE(asin);
    VECTORIZE(atan);
    VECTORIZE2(atan2);
    VECTORIZE(cosh);
    VECTORIZE(sinh);
    VECTORIZE(tanh);
    VECTORIZE(acosh);
    VECTORIZE(asinh);
    VECTORIZE(atanh);
    VECTORIZE(exp);
    VECTORIZE(log);
    VECTORIZE(log2);
    VECTORIZE(log10);
    VECTORIZE(erf);
    VECTORIZE(erfc);
    VECTORIZE(tgamma);
    VECTORIZE(ceil);
    VECTORIZE(floor);
    VECTORIZE(round);
    VECTORIZE(abs);
    VECTORIZE(clamp);
    VECTORIZE_REN(bessel_j0, j0);
    VECTORIZE_REN(bessel_j1, j1);
    VECTORIZE_REN(bessel_y0, y0);
    VECTORIZE_REN(bessel_y1, y1);

    #ifndef NO_GSL
    VECTORIZE_REN(bessel_i0, gsl_sf_bessel_I0);
    VECTORIZE_REN(bessel_i1, gsl_sf_bessel_I1);
    VECTORIZE_REN(bessel_k0, gsl_sf_bessel_K0);
    VECTORIZE_REN(bessel_k1, gsl_sf_bessel_K1);
    #endif

    #undef VECTORIZE

    // Create a range of n steps from i to j (inclusive)
    template<typename T, typename U = T>
    vec1d rgen(T i, U j, uint_t n) {
        if (n == 1) {
            return {i};
        } else {
            vec1d v(n);
            double dx = (j-i)/double(n-1);
            for (uint_t k : range(uint_t(n))) {
                v.safe[k] = i + k*dx;
            }

            return v;
        }
    }

    // Create a range of n logarithmic steps from i to j (inclusive)
    template<typename T, typename U = T>
    vec1d rgen_log(T i, U j, uint_t n) {
        phypp_check(i > 0 && j > 0, "'rgen_log(a,b,n)' needs a strictly positive value for 'a' and 'b' "
            "(got ", i, " and ", j, ")");

        if (n == 1) {
            return {i};
        } else {
            vec1d v(n);
            double dx = log10(j/i)/double(n-1);
            for (uint_t k : range(uint_t(n))) {
                v.safe[k] = i*e10(k*dx);
            }

            return v;
        }
    }

    // Create a range from i to j (inclusive) with step s
    // If the range is not exactly dividable by s, the last value is forced to equal j.
    // The range will always contain both i and j.
    template<typename T, typename U, typename V>
    vec1d rgen_step(T i, U j, V s) {
        vec1d ret;

        phypp_check(s > 0 && is_finite(s),
            "'rgen_step(a,b,s)' needs a strictly positive and finite value for 's' (got ", s, ")");

        uint_t n;
        double step = s;
        if (i < T(j)) {
            n = round((j - i)/s) + 1;
        } else {
            n = round((i - j)/s) + 1;
            step *= -1.0;
        }

        ret = dindgen(n)*step + i;

        if (n > 1) {
            ret.back() = j;
        } else if (i - j != 0.0) {
            ret.push_back(j);
        }

        return ret;
    }
}

#endif
