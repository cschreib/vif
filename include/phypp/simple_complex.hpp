#ifndef SIMPLE_COMPLEX_HPP
#define SIMPLE_COMPLEX_HPP

#ifdef USE_PHYPP_COMPLEX

namespace phypp {
    template<typename T>
    struct complex {
        T re, im;

        template<typename U>
        void operator *= (U&& c) {
            *this = (*this)*std::forward<U>(c);
        }

        template<typename U>
        void operator /= (U&& c) {
            *this = (*this)/std::forward<U>(c);
        }

        template<typename U>
        void operator += (U&& c) {
            *this = (*this) + std::forward<U>(c);
        }

        template<typename U>
        void operator -= (U&& c) {
            *this = (*this) - std::forward<U>(c);
        }
    };

    template<typename T>
    T norm_sqr(complex<T> c) {
        return c.re*c.re + c.im*c.im;
    }

    template<typename T>
    T norm(complex<T> c) {
        return sqrt(norm_sqr(c));
    }

    template<typename T, typename U>
    auto operator* (complex<T> c1, complex<U> c2) -> complex<decltype(c1.re*c2.re)> {
        return {c1.re*c2.re - c1.im*c2.im, c1.re*c2.im + c1.im*c2.re};
    }

    template<typename T, typename U, typename enable =
        typename std::enable_if<std::is_arithmetic<U>::value>::type>
    auto operator* (complex<T> c, U a) -> complex<decltype(c.re*a)> {
        return {c.re*a, c.im*a};
    }

    template<typename T, typename U, typename enable =
        typename std::enable_if<std::is_arithmetic<U>::value>::type>
    auto operator* (U a, complex<T> c) -> complex<decltype(a*c.re)> {
        return {a*c.re, a*c.im};
    }

    template<typename T, typename U>
    auto operator/ (complex<T> c1, complex<U> c2) -> complex<decltype(c1.re*c2.re)> {
        auto tmp = 1.0/norm_sqr(c2);
        return {c1.re*(c2.re*tmp) + c1.im*(c2.im*tmp), -c1.re*(c2.im*tmp) + c1.im*(c2.re*tmp)};
    }

    template<typename T, typename U, typename enable =
        typename std::enable_if<std::is_arithmetic<U>::value>::type>
    auto operator/ (complex<T> c, U a) -> complex<decltype(c.re/a)> {
        return {c.re/a, c.im/a};
    }

    template<typename T, typename U, typename enable =
        typename std::enable_if<std::is_arithmetic<U>::value>::type>
    auto operator/ (U a, complex<T> c) -> complex<decltype(a/c.re)> {
        auto tmp = 1.0/norm_sqr(c);
        return {a*(c.re*tmp), -a*(c.im*tmp)};
    }

    template<typename T, typename U>
    auto operator+ (complex<T> c1, complex<U> c2) -> complex<decltype(c1.re + c2.re)> {
        return {c1.re + c2.re, c1.im + c2.im};
    }

    template<typename T, typename U, typename enable =
        typename std::enable_if<std::is_arithmetic<U>::value>::type>
    auto operator+ (complex<T> c, U a) -> complex<decltype(c.re + a)> {
        return {c.re + a, c.im};
    }

    template<typename T, typename U, typename enable =
        typename std::enable_if<std::is_arithmetic<U>::value>::type>
    auto operator+ (U a, complex<T> c) -> complex<decltype(a + c.re)> {
        return {c.re + a, c.im};
    }

    template<typename T, typename U>
    auto operator- (complex<T> c1, complex<U> c2) -> complex<decltype(c1.re - c2.re)> {
        return {c1.re - c2.re, c1.im - c2.im};
    }

    template<typename T, typename U, typename enable =
        typename std::enable_if<std::is_arithmetic<U>::value>::type>
    auto operator- (complex<T> c, U a) -> complex<decltype(c.re - a)> {
        return {c.re - a, c.im};
    }

    template<typename T, typename U, typename enable =
        typename std::enable_if<std::is_arithmetic<U>::value>::type>
    auto operator- (U a, complex<T> c) -> complex<decltype(a - c.re)> {
        return {a - c.re, -c.im};
    }

    template<typename T>
    complex<T> operator+ (complex<T> c) {
        return c;
    }

    template<typename T>
    complex<T> operator- (complex<T> c) {
        return {-c.re, -c.im};
    }
}

#define MAKE_TYPEDEFS(N) \
    using vec##N##cf  = vec<N, phypp::complex<float>>; \
    using vec##N##cd  = vec<N, phypp::complex<double>>;
#else

#include <complex>

#define MAKE_TYPEDEFS(N) \
    using vec##N##cf  = vec<N, std::complex<float>>; \
    using vec##N##cd  = vec<N, std::complex<double>>;

#endif

MAKE_TYPEDEFS(1)
MAKE_TYPEDEFS(2)
MAKE_TYPEDEFS(3)
MAKE_TYPEDEFS(4)
MAKE_TYPEDEFS(5)
MAKE_TYPEDEFS(6)

#undef MAKE_TYPEDEFS

#endif
