#ifndef PHYPP_MATH_INTERPOLATE_HPP
#define PHYPP_MATH_INTERPOLATE_HPP

#include "phypp/core/vec.hpp"
#include "phypp/core/range.hpp"
#include "phypp/core/error.hpp"
#include "phypp/utility/generic.hpp"
#include "phypp/math/base.hpp"

namespace phypp {
    inline double interpolate(double y1, double y2, double x1, double x2, double x) {
        double a = (x - x1)/(x2 - x1);
        return y1 + (y2 - y1)*a;
    }

    inline std::pair<double,double> interpolate_err(double y1, double y2, double e1, double e2, double x1, double x2, double x) {
        double a = (x - x1)/(x2 - x1);
        return std::make_pair(y1 + (y2 - y1)*a, sqrt(sqr(e1*(1-a)) + sqr(e2*a)));
    }

    // Perform linear interpolation of data 'y' of position 'x' at new positions 'nx'.
    // Assumes that the arrays only contain finite elements, and that 'x' is properly sorted. If one of
    // the arrays contains special values (NaN, inf, ...), all the points that would use these values
    // will be contaminated. If 'x' is not properly sorted, the result will simply be wrong.
    template<std::size_t DI = 1, std::size_t DX = 1, typename TypeX2 = double, typename TypeY = double,
        typename TypeX1 = double>
    auto interpolate(const vec<DI,TypeY>& y, const vec<DI,TypeX1>& x, const vec<DX,TypeX2>& nx) ->
        vec<DX,decltype(y[0]*x[0])> {

        using rtypey = meta::rtype_t<TypeY>;
        using rtypex = meta::rtype_t<TypeX1>;

        phypp_check(y.size() == x.size(),
            "interpolate: 'x' and 'y' arrays must contain the same number of elements");
        phypp_check(y.size() >= 2,
            "interpolate: 'x' and 'y' arrays must contain at least 2 elements");

        uint_t nmax = x.size();
        vec<DX,decltype(y[0]*x[0])> r; r.reserve(nx.size());
        for (auto& tx : nx) {
            uint_t low = lower_bound(tx, x);

            rtypey ylow, yup;
            rtypex xlow, xup;
            if (low != npos) {
                if (low != nmax-1) {
                    ylow = y.safe[low]; yup = y.safe[low+1];
                    xlow = x.safe[low]; xup = x.safe[low+1];
                } else {
                    ylow = y.safe[low-1]; yup = y.safe[low];
                    xlow = x.safe[low-1]; xup = x.safe[low];
                }
            } else {
                ylow = y.safe[0]; yup = y.safe[1];
                xlow = x.safe[0]; xup = x.safe[1];
            }

            r.push_back(ylow + (yup - ylow)*(tx - xlow)/(xup - xlow));
        }

        return r;
    }

    // Perform linear interpolation of data 'y' of position 'x' at new position 'nx'.
    // Assumes that the arrays only contain finite elements, and that 'x' is properly sorted.
    // If one of the arrays contains special values (NaN, inf, ...), all the points that would use
    // these values will be contaminated. If 'x' is not properly sorted, the result will simply be
    // wrong.
    template<std::size_t DI = 1, typename TypeY = double, typename TypeX = double,
        typename T = double,
        typename enable = typename std::enable_if<!meta::is_vec<T>::value>::type>
    auto interpolate(const vec<DI,TypeY>& y, const vec<DI,TypeX>& x, const T& nx) ->
        decltype(y[0]*x[0]) {

        using rtypey = meta::rtype_t<TypeY>;
        using rtypex = meta::rtype_t<TypeX>;

        phypp_check(y.dims == x.dims, "incompatible dimensions between X and Y arrays (", x.dims,
            " vs. ", y.dims, ")");
        phypp_check(y.size() >= 2, "interpolated data must contain at least 2 elements");

        uint_t nmax = x.size();
        uint_t low = lower_bound(nx, x);

        rtypey ylow, yup;
        rtypex xlow, xup;
        if (low != npos) {
            if (low != nmax-1) {
                ylow = y.safe[low]; yup = y.safe[low+1];
                xlow = x.safe[low]; xup = x.safe[low+1];
            } else {
                ylow = y.safe[low-1]; yup = y.safe[low];
                xlow = x.safe[low-1]; xup = x.safe[low];
            }
        } else {
            ylow = y.safe[0]; yup = y.safe[1];
            xlow = x.safe[0]; xup = x.safe[1];
        }

        return ylow + (yup - ylow)*(nx - xlow)/(xup - xlow);
    }

    // Perform linear interpolation of data 'y' of position 'x' at new positions 'nx'.
    // Assumes that the arrays only contain finite elements, and that 'x' is properly sorted. If one of
    // the arrays contains special values (NaN, inf, ...), all the points that would use these values
    // will be contaminated. If 'x' is not properly sorted, the result will simply be wrong.
    template<std::size_t DI = 1, std::size_t DX = 1, typename TypeX2 = double, typename TypeY = double,
        typename TypeE = double, typename TypeX1 = double>
    auto interpolate(const vec<DI,TypeY>& y, const vec<DI,TypeE>& e, const vec<DI,TypeX1>& x,
        const vec<DX,TypeX2>& nx) ->
        std::pair<vec<DX,decltype(y[0]*x[0])>,vec<DX,decltype(y[0]*x[0])>> {

        using rtypey = meta::rtype_t<TypeY>;
        using rtypee = meta::rtype_t<TypeE>;
        using rtypex = meta::rtype_t<TypeX1>;

        phypp_check(y.dims == x.dims, "incompatible dimensions between X and Y arrays (", x.dims,
            " vs. ", y.dims, ")");
        phypp_check(y.size() >= 2, "interpolated data must contain at least 2 elements");

        uint_t nmax = x.size();
        std::pair<vec<DX,decltype(y[0]*x[0])>,vec<DX,decltype(y[0]*x[0])>> p;
        p.first.reserve(nx.size());
        p.second.reserve(nx.size());
        for (auto& tx : nx) {
            uint_t low = lower_bound(tx, x);

            rtypey ylow, yup;
            rtypee elow, eup;
            rtypex xlow, xup;
            if (low != npos) {
                if (low != nmax-1) {
                    ylow = y.safe[low]; yup = y.safe[low+1];
                    elow = e.safe[low]; eup = e.safe[low+1];
                    xlow = x.safe[low]; xup = x.safe[low+1];
                } else {
                    ylow = y.safe[low-1]; yup = y.safe[low];
                    elow = e.safe[low-1]; eup = e.safe[low];
                    xlow = x.safe[low-1]; xup = x.safe[low];
                }
            } else {
                ylow = y.safe[0]; yup = y.safe[1];
                elow = e.safe[0]; eup = e.safe[1];
                xlow = x.safe[0]; xup = x.safe[1];
            }

            double a = (tx - xlow)/(xup - xlow);
            p.first.push_back(ylow + (yup - ylow)*a);
            p.second.push_back(sqrt(sqr(elow*(1-a)) + sqr(eup*a)));
        }

        return p;
    }

    // Raw bilinear interpolation (or extrapolation) between four values
    inline double bilinear(double v1, double v2, double v3, double v4, double x, double y) {
        return v1*(1.0-x)*(1.0-y) + v2*(1.0-x)*y + v3*x*(1.0-y) + v4*x*y;
    }

    template<typename Type>
    meta::rtype_t<Type> bilinear(const vec<2,Type>& map, double x, double y) {
        int_t tix = floor(x);
        int_t tiy = floor(y);
        double dx = x - tix;
        double dy = y - tiy;

        if (tix < 0) {
            tix = 0;
            dx = x - tix;
        } else if (uint_t(tix) >= map.dims[0]-1) {
            tix = map.dims[0]-2;
            dx = x - tix;
        }

        if (tiy < 0) {
            tiy = 0;
            dy = y - tiy;
        } else if (uint_t(tiy) >= map.dims[1]-1) {
            tiy = map.dims[1]-2;
            dy = y - tiy;
        }

        uint_t ix = tix;
        uint_t iy = tiy;

        return bilinear(map.safe(ix,iy), map.safe(ix,iy+1),
                        map.safe(ix+1,iy), map.safe(ix+1,iy+1), dx, dy);
    }

    template<typename Type>
    meta::rtype_t<Type> bilinear(const vec<2,Type>& map, const vec1d& mx,
        const vec1d& my, double x, double y) {

        phypp_check(map.dims[0] == mx.size(), "incompatible size of MAP and MX (", map.dims,
            " vs. ", mx.size(), ")");
        phypp_check(map.dims[1] == my.size(), "incompatible size of MAP and MY (", map.dims,
            " vs. ", my.size(), ")");

        double ux = interpolate(dindgen(mx.size()), mx, x);
        double uy = interpolate(dindgen(my.size()), my, y);
        return bilinear(map, ux, uy);
    }

    template<typename Type, std::size_t D, typename TX, typename TY>
    vec<D,meta::rtype_t<Type>> bilinear(const vec<2,Type>& map, const vec1d& mx,
        const vec1d& my, const vec<D,TX>& x, const vec<D,TY>& y) {

        phypp_check(map.dims[0] == mx.size(), "incompatible size of MAP and MX (", map.dims,
            " vs. ", mx.size(), ")");
        phypp_check(map.dims[1] == my.size(), "incompatible size of MAP and MY (", map.dims,
            " vs. ", my.size(), ")");
        phypp_check(x.dims == y.dims, "incompatible dimensions for X and Y (", x.dims,
            " vs. ", y.size(), ")");

        vec<D,meta::rtype_t<Type>> v(x.dims);

        auto ux = interpolate(dindgen(mx.size()), mx, x);
        auto uy = interpolate(dindgen(my.size()), my, y);

        for (uint_t i : range(x)) {
            v.safe[i] = bilinear(map, ux.safe[i], uy.safe[i]);
        }

        return v;
    }


    template<typename Type, typename TypeD = double>
    meta::rtype_t<Type> bilinear_strict(const vec<2,Type>& map, double x, double y,
        TypeD def = 0.0) {

        int_t tix = floor(x);
        int_t tiy = floor(y);

        if (tix < 0 || uint_t(tix) >= map.dims[0]-1 || tiy < 0 || uint_t(tiy) >= map.dims[1]-1) {
            return def;
        }

        double dx = x - tix;
        double dy = y - tiy;
        uint_t ix = tix;
        uint_t iy = tiy;

        return bilinear(map.safe(ix,iy), map.safe(ix,iy+1),
                        map.safe(ix+1,iy), map.safe(ix+1,iy+1), dx, dy);
    }

    template<typename Type>
    vec<2,meta::rtype_t<Type>> rebin(const vec<2,Type>& map, const vec1d& mx,
        const vec1d& my, const vec1d& x, const vec1d& y) {

        phypp_check(map.dims[0] == mx.size(), "incompatible size of MAP and MX (", map.dims,
            " vs. ", mx.size(), ")");
        phypp_check(map.dims[1] == my.size(), "incompatible size of MAP and MY (", map.dims,
            " vs. ", my.size(), ")");

        vec<2,meta::rtype_t<Type>> v(x.size(), y.size());

        vec1d ux = interpolate(dindgen(mx.size()), mx, x);
        vec1d uy = interpolate(dindgen(my.size()), my, y);

        for (uint_t ix : range(ux))
        for (uint_t iy : range(uy)) {
            v.safe(ix,iy) = bilinear(map, ux.safe[ix], uy.safe[iy]);
        }

        return v;
    }

    // Perform cubic interpolation of data 'y' of position 'x' at new position 'nx'.
    // Assumes that the arrays only contain finite elements, and that 'x' is properly sorted.
    // If one of the arrays contains special values (NaN, inf, ...), all the points that would use
    // these values will be contaminated. If 'x' is not properly sorted, the result will simply be
    // wrong.
    template <std::size_t D, typename TN, typename TX, typename TY>
    vec<D,meta::rtype_t<TY>> interpolate_3spline(const vec<1,TY>& y, const vec<1,TX>& x,
        const vec<D,TN>& xn) {

        vec<D,meta::rtype_t<TY>> yn(xn.dims);

        phypp_check(y.size() == x.size(),
            "'x' and 'y' arrays must contain the same number of elements");
        phypp_check(y.size() >= 2,
            "'x' and 'y' arrays must contain at least 2 elements");

        const uint_t n = y.size()-1;
        vec1d b(n+1), d(n+1);
        vec1d h(n);
        for (uint_t i : range(h)) {
            h[i] = x[i+1]-x[i];
        }

        vec1d al(n);
        for (uint_t i : range(1, n)) {
            al[i] = 3.0*(y[i+1] - y[i])/h[i] - 3.0*(y[i] - y[i-1])/h[i];
        }

        vec1d c(n+1), l(n+1), mu(n+1), z(n+1);
        l[0] = 1;
        for (uint_t i : range(1, n)) {
            l[i] = 2.0*(x[i+1] - x[i-1]) - h[i-1]*mu[i-1];
            mu[i] = h[i]/l[i];
            z[i] = (al[i] - h[i-1]*z[i-1])/l[i];
        }

        l[n] = 1; {
            uint_t i = n;
            do {
                --i;
                c[i] = z[i] - mu[i]*c[i+1];
                b[i] = (y[i+1]-y[i])/h[i] - (1.0/3.0)*h[i]*(c[i+1] + 2*c[i]);
                d[i] = (1.0/3.0)*(c[i+1] - c[i])/h[i];
            } while (i != 0);
        }

        b[n] = b[n-1] + 2*c[n-1]*(x[n]-x[n-1]) + 3*d[n-1]*sqr(x[n]-x[n-1]);

        for (uint_t i : range(xn)) {
            uint_t k = lower_bound(xn.safe[i], x);
            if (k == npos) {
                double th = xn[i] - x[0];
                yn[i] = y[0] + b[0]*th;
            } else {
                double th = xn[i] - x[k];
                yn[i] = y[k] + b[k]*th + c[k]*th*th + d[k]*th*th*th;
            }
        }

        return yn;
    }

    // Perform cubic interpolation of the regularly gridded data 'y' at the new positions 'nx'.
    // Assumes that the data array only contains finite elements. If this is not the case,
    // all the points that would use these values will be contaminated.
    template <std::size_t D, typename TN, typename TY>
    vec<D,meta::rtype_t<TY>> interpolate_3spline(const vec<1,TY>& y, const vec<D,TN>& xn) {
        vec<D,meta::rtype_t<TY>> yn(xn.dims);

        phypp_check(y.size() >= 2, "'y' array must contain at least 2 elements");

        const uint_t n = y.size()-1;
        vec1d b(n+1), d(n+1);

        vec1d al(n);
        for (uint_t i : range(1, n)) {
            al[i] = 3.0*(y[i+1] - y[i]) - 3.0*(y[i] - y[i-1]);
        }

        vec1d c(n+1), l(n+1), mu(n+1), z(n+1);
        l[0] = 1;
        for (uint_t i : range(1, n)) {
            l[i] = 4.0 - mu[i-1];
            mu[i] = 1.0/l[i];
            z[i] = (al[i] - z[i-1])/l[i];
        }

        l[n] = 1; {
            uint_t i = n;
            do {
                --i;
                c[i] = z[i] - mu[i]*c[i+1];
                b[i] = (y[i+1]-y[i]) - (1.0/3.0)*(c[i+1] + 2*c[i]);
                d[i] = (1.0/3.0)*(c[i+1] - c[i]);
            } while (i != 0);
        }

        b[n] = b[n-1] + 2*c[n-1] + 3*d[n-1];

        for (uint_t i : range(xn)) {
            int_t k = floor(xn.safe[i]);
            if (k < 0) {
                double th = xn[i];
                yn[i] = y[0] + b[0]*th;
            } else {
                if (k > n) k = n;
                double th = double(xn[i]) - k;
                yn[i] = y[k] + b[k]*th + c[k]*th*th + d[k]*th*th*th;
            }
        }

        return yn;
    }

    // Perform cubic interpolation of the regularly gridded data 'y' at the new positions 'nx'.
    // Assumes that the data array only contains finite elements. If this is not the case,
    // all the points that would use these values will be contaminated.
    template <typename TN, typename TY>
    meta::rtype_t<TY> interpolate_3spline(const vec<1,TY>& y, const TN& xn) {
        phypp_check(y.size() >= 2, "'y' array must contain at least 2 elements");

        const uint_t n = y.size()-1;
        vec1d b(n+1), d(n+1);

        vec1d al(n);
        for (uint_t i : range(1, n)) {
            al[i] = 3.0*(y[i+1] - y[i]) - 3.0*(y[i] - y[i-1]);
        }

        vec1d c(n+1), l(n+1), mu(n+1), z(n+1);
        l[0] = 1;
        for (uint_t i : range(1, n)) {
            l[i] = 4.0 - mu[i-1];
            mu[i] = 1.0/l[i];
            z[i] = (al[i] - z[i-1])/l[i];
        }

        l[n] = 1; {
            uint_t i = n;
            do {
                --i;
                c[i] = z[i] - mu[i]*c[i+1];
                b[i] = (y[i+1]-y[i]) - (1.0/3.0)*(c[i+1] + 2*c[i]);
                d[i] = (1.0/3.0)*(c[i+1] - c[i]);
            } while (i != 0);
        }

        b[n] = b[n-1] + 2*c[n-1] + 3*d[n-1];

        int_t k = floor(xn);
        if (k < 0) {
            double th = xn;
            return y[0] + b[0]*th;
        } else {
            if (k > n) k = n;
            double th = double(xn) - k;
            return y[k] + b[k]*th + c[k]*th*th + d[k]*th*th*th;
        }
    }

    // Perform bicubic interpolation of the regularly gridded data 'map' at the new positions
    // 'x' and 'y'.
    // Assumes that the data array only contains finite elements. If this is not the case,
    // all the points that would use these values will be contaminated.
    template<typename Type, typename TypeD = double>
    meta::rtype_t<Type> bicubic_strict(const vec<2,Type>& map, double x, double y,
        TypeD def = 0.0) {

        int_t tix = floor(x);
        int_t tiy = floor(y);

        if (tix <= 0 || uint_t(tix) >= map.dims[0]-2 || tiy <= 0 || uint_t(tiy) >= map.dims[1]-2) {
            return def;
        }

        vec<1,meta::rtype_t<Type>> tmpx(4);
        vec<1,meta::rtype_t<Type>> tmpy(4);
        for (int_t ix : range(tmpx)) {
            uint_t uix = ix + tix - 1;
            for (uint_t iy : range(tmpy)) {
                uint_t uiy = iy + tiy - 1;
                tmpy.safe[iy] = map.safe(uix, uiy);
            }

            tmpx.safe[ix] = interpolate_3spline(tmpy, y - tiy + 1);
        }

        return interpolate_3spline(tmpx, x - tix + 1);
    }
}

#endif
