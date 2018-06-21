#ifndef PHYPP_MATH_MATH_HPP
#define PHYPP_MATH_MATH_HPP

#include <cmath>
#include <limits>
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
    static constexpr const double dpi = 3.14159265358979323;
    static constexpr const float  fpi = 3.14159265358979323;
    static constexpr const double ln10 = 2.30258509299404568;

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
    auto e10(const T& t) -> decltype(exp(ln10*t)) {
        return exp(ln10*t);
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

    PHYPP_VECTORIZE(sqrt);
    PHYPP_VECTORIZE(sqr);
    PHYPP_VECTORIZE(invsqr);
    PHYPP_VECTORIZE(pow);
    PHYPP_VECTORIZE(fmod);
    PHYPP_VECTORIZE(cos);
    PHYPP_VECTORIZE(sin);
    PHYPP_VECTORIZE(tan);
    PHYPP_VECTORIZE(acos);
    PHYPP_VECTORIZE(asin);
    PHYPP_VECTORIZE(atan);
    PHYPP_VECTORIZE2(atan2);
    PHYPP_VECTORIZE(cosh);
    PHYPP_VECTORIZE(sinh);
    PHYPP_VECTORIZE(tanh);
    PHYPP_VECTORIZE(acosh);
    PHYPP_VECTORIZE(asinh);
    PHYPP_VECTORIZE(atanh);
    PHYPP_VECTORIZE(exp);
    PHYPP_VECTORIZE(log);
    PHYPP_VECTORIZE(log2);
    PHYPP_VECTORIZE(log10);
    PHYPP_VECTORIZE(erf);
    PHYPP_VECTORIZE(erfc);
    PHYPP_VECTORIZE(tgamma);
    PHYPP_VECTORIZE(ceil);
    PHYPP_VECTORIZE(floor);
    PHYPP_VECTORIZE(round);
    PHYPP_VECTORIZE(abs);
    PHYPP_VECTORIZE(clamp);
    PHYPP_VECTORIZE_REN(bessel_j0, j0);
    PHYPP_VECTORIZE_REN(bessel_j1, j1);
    PHYPP_VECTORIZE_REN(bessel_y0, y0);
    PHYPP_VECTORIZE_REN(bessel_y1, y1);

    #ifndef NO_GSL
    PHYPP_VECTORIZE_REN(bessel_i0, gsl_sf_bessel_I0);
    PHYPP_VECTORIZE_REN(bessel_i1, gsl_sf_bessel_I1);
    PHYPP_VECTORIZE_REN(bessel_k0, gsl_sf_bessel_K0);
    PHYPP_VECTORIZE_REN(bessel_k1, gsl_sf_bessel_K1);
    #endif

    // Create a range of n steps from i to j (inclusive)
    template<typename T, typename U = T>
    vec1d rgen(T i, U j, uint_t n) {
        if (n == 1) {
            return vec1d({static_cast<double>(i)});
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

        using ctype = decltype(i*j*s);

        uint_t n;
        double step = s;
        if (ctype(i) < ctype(j)) {
            n = round((ctype(j) - ctype(i))/s) + 1;
        } else {
            n = round((ctype(i) - ctype(j))/s) + 1;
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
