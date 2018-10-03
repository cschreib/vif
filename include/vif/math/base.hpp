#ifndef VIF_MATH_MATH_HPP
#define VIF_MATH_MATH_HPP

#include <cmath>
#include <limits>
#ifndef NO_GSL
#include <gsl/gsl_sf_bessel.h>
#endif
#include "vif/core/vec.hpp"
#include "vif/core/range.hpp"
#include "vif/core/error.hpp"

namespace vif {
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

    float fast_exp(float x) {
        // Implementation based on:
        // https://stackoverflow.com/a/10792321/1565581

        static_assert(sizeof(float) == sizeof(std::int32_t), "size of float and int32_t don't match");

        if (x < -87.0f || x > 88.0f) return 0.0f;

        // exp(x) = 2^i * 2^f; i = floor (log2(e) * x), 0 <= f <= 1
        float t = x * 1.442695041f;
        float fi = floor(t);
        float f = t - fi;
        std::int32_t i = (std::int32_t)fi;

        // compute 2^f
        f = (0.3371894346f*f + 0.657636276f)*f + 1.00172476f;

        // scale by 2^i
        std::int32_t k;
        std::memcpy(&k, &f, sizeof(float));
        k += (i << 23);
        std::memcpy(&f, &k, sizeof(float));

        return f;
    }

    VIF_VECTORIZE(sqrt);
    VIF_VECTORIZE(sqr);
    VIF_VECTORIZE(invsqr);
    VIF_VECTORIZE(pow);
    VIF_VECTORIZE(fmod);
    VIF_VECTORIZE(cos);
    VIF_VECTORIZE(sin);
    VIF_VECTORIZE(tan);
    VIF_VECTORIZE(acos);
    VIF_VECTORIZE(asin);
    VIF_VECTORIZE(atan);
    VIF_VECTORIZE2(atan2);
    VIF_VECTORIZE(cosh);
    VIF_VECTORIZE(sinh);
    VIF_VECTORIZE(tanh);
    VIF_VECTORIZE(acosh);
    VIF_VECTORIZE(asinh);
    VIF_VECTORIZE(atanh);
    VIF_VECTORIZE(exp);
    VIF_VECTORIZE(fast_exp);
    VIF_VECTORIZE(log);
    VIF_VECTORIZE(log2);
    VIF_VECTORIZE(log10);
    VIF_VECTORIZE(erf);
    VIF_VECTORIZE(erfc);
    VIF_VECTORIZE(tgamma);
    VIF_VECTORIZE(ceil);
    VIF_VECTORIZE(floor);
    VIF_VECTORIZE(round);
    VIF_VECTORIZE(abs);
    VIF_VECTORIZE(clamp);
    VIF_VECTORIZE_REN(bessel_j0, j0);
    VIF_VECTORIZE_REN(bessel_j1, j1);
    VIF_VECTORIZE_REN(bessel_y0, y0);
    VIF_VECTORIZE_REN(bessel_y1, y1);

    #ifndef NO_GSL
    VIF_VECTORIZE_REN(bessel_i0, gsl_sf_bessel_I0);
    VIF_VECTORIZE_REN(bessel_i1, gsl_sf_bessel_I1);
    VIF_VECTORIZE_REN(bessel_k0, gsl_sf_bessel_K0);
    VIF_VECTORIZE_REN(bessel_k1, gsl_sf_bessel_K1);
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
        vif_check(i > 0 && j > 0, "'rgen_log(a,b,n)' needs a strictly positive value for 'a' and 'b' "
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

        vif_check(s > 0 && is_finite(s),
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

        ret = indgen<double>(n)*step + i;

        if (n > 1) {
            ret.back() = j;
        } else if (i - j != 0.0) {
            ret.push_back(j);
        }

        return ret;
    }
}

#endif
