#ifndef PHYPP_MATH_HISTOGRAM_HPP
#define PHYPP_MATH_HISTOGRAM_HPP

#include "phypp/core/vec.hpp"
#include "phypp/core/range.hpp"
#include "phypp/core/error.hpp"
#include "phypp/utility/generic.hpp"
#include "phypp/math/base.hpp"

namespace phypp {
    template<typename T>
    vec<2,T> make_bins(T mi, T ma) {
        return {{mi}, {ma}};
    }

    template<typename T>
    vec2d make_bins(T mi, T ma, uint_t n) {
        vec2d b(2,n);
        double d = (ma - mi)/double(n);
        for (uint_t i : range(n)) {
            b.safe(0,i) = mi + i*d;
            b.safe(1,i) = mi + (i+1)*d;
        }

        return b;
    }

    template<typename T = double>
    vec<2,T> make_bins(const vec<1,T>& v) {
        vec<2,T> b(2, v.size()-1);
        for (uint_t i : range(b.dims[1])) {
            b.safe(0,i) = v.safe[i];
            b.safe(1,i) = v.safe[i+1];
        }

        return b;
    }

    template<typename T, typename B = double>
    bool in_bin(T t, const vec<1,B>& b) {
        phypp_check(b.dims[0] == 2, "can only be called on a single bin "
            "vector (expected dims=[2], got dims=[", b.dims, "])");

        return t >= b.safe[0] && t < b.safe[1];
    }

    template<std::size_t Dim, typename Type, typename B = double>
    vec<Dim,bool> in_bin(const vec<Dim,Type>& v, const vec<1,B>& b) {
        phypp_check(b.dims[0] == 2, "can only be called on a single bin "
            "vector (expected dims=[2], got dims=[", b.dims, "])");

        auto low = b.safe[0];
        auto up  = b.safe[1];
        vec<Dim,bool> res(v.dims);
        for (uint_t i : range(v)) {
            res.safe[i] = v.safe[i] >= low && v.safe[i] < up;
        }

        return res;
    }

    template<typename T, typename B = double>
    bool in_bin(T t, const vec<2,B>& b, uint_t ib) {
        phypp_check(b.dims[0] == 2, "can only be called on a single bin "
            "vector (expected dims=[2], got dims=[", b.dims, "])");
        phypp_check(ib < b.dims[1], "bin index is out of bounds ",
            "(", ib, " vs. ", b.dims[1], ")");

        return t >= b.safe(0,ib) && t < b.safe(1,ib);
    }

    template<std::size_t Dim, typename Type, typename B = double>
    vec<Dim,bool> in_bin(const vec<Dim,Type>& v, const vec<2,B>& b, uint_t ib) {
        phypp_check(b.dims[0] == 2, "B is not a bin vector "
            "(expected dims=[2, ...], got dims=[", b.dims, "])");
        phypp_check(ib < b.dims[1], "bin index is out of bounds ",
            "(", ib, " vs. ", b.dims[1], ")");

        auto low = b.safe(0,ib);
        auto up  = b.safe(1,ib);
        vec<Dim,bool> res(v.dims);
        for (uint_t i : range(v)) {
            res.safe[i] = v.safe[i] >= low && v.safe[i] < up;
        }

        return res;
    }

    template<typename T, typename B = double>
    bool in_bin_open(T t, const vec<2,B>& b, uint_t ib) {
        phypp_check(b.dims[0] == 2, "can only be called on a single bin "
            "vector (expected dims=[2], got dims=[", b.dims, "])");
        phypp_check(ib < b.dims[1], "bin index is out of bounds ",
            "(", ib, " vs. ", b.dims[1], ")");

        if (ib == 0) {
            return t < b.safe(1,ib);
        } else if (ib == b.dims[1]-1) {
            return t >= b.safe(0,ib);
        } else {
            return t >= b.safe(0,ib) && t < b.safe(1,ib);
        }
    }

    template<std::size_t Dim, typename Type, typename B = double>
    vec<Dim,bool> in_bin_open(const vec<Dim,Type>& v, const vec<2,B>& b, uint_t ib) {
        phypp_check(b.dims[0] == 2, "B is not a bin vector "
            "(expected dims=[2, ...], got dims=[", b.dims, "])");
        phypp_check(ib < b.dims[1], "bin index is out of bounds ",
            "(", ib, " vs. ", b.dims[1], ")");

        auto low = b.safe(0,ib);
        auto up  = b.safe(1,ib);
        vec<Dim,bool> res(v.dims);
        if (ib == 0) {
            for (uint_t i : range(v)) {
                res.safe[i] = v.safe[i] < up;
            }
        } else if (ib == b.dims[1]-1) {
            for (uint_t i : range(v)) {
                res.safe[i] = v.safe[i] >= low;
            }
        } else {
            for (uint_t i : range(v)) {
                res.safe[i] = v.safe[i] >= low && v.safe[i] < up;
            }
        }

        return res;
    }

    template<typename Type>
    auto bin_center(const vec<2,Type>& b) -> vec<1,decltype(0.5*b[0])> {
        phypp_check(b.dims[0] == 2, "B is not a bin vector "
            "(expected dims=[2, ...], got dims=[", b.dims, "])");

        return 0.5*(b.safe(1,_) + b.safe(0,_));
    }

    template<typename Type>
    auto bin_center(const vec<1,Type>& b) -> decltype(0.5*b[0]) {
        phypp_check(b.dims[0] == 2, "B is not a bin vector "
            "(expected dims=2, got dims=", b.dims, ")");

        return 0.5*(b.safe[1] + b.safe[0]);
    }

    template<typename Type>
    vec<1,meta::rtype_t<Type>> bin_width(const vec<2,Type>& b) {
        phypp_check(b.dims[0] == 2, "B is not a bin vector "
            "(expected dims=[2, ...], got dims=[", b.dims, "])");

        return b.safe(1,_) - b.safe(0,_);
    }

    template<typename Type>
    meta::rtype_t<Type> bin_width(const vec<1,Type>& b) {
        phypp_check(b.dims[0] == 2, "B is not a bin vector "
            "(expected dims=2, got dims=", b.dims, ")");

        return b.safe[1] - b.safe[0];
    }

    template<std::size_t Dim, typename Type, typename TypeB>
    vec1u histogram(const vec<Dim,Type>& data, const vec<2,TypeB>& bins) {
        phypp_check(bins.dims[0] == 2, "can only be called with a bin vector (expected "
            "dims=[2,...], got dims=[", bins.dims, "])");

        using rtype = meta::rtype_t<Type>;
        vec<Dim,rtype> tmp = data;

        uint_t nbin = bins.dims[1];
        vec1u counts(nbin);

        auto first = tmp.data.begin();
        for (uint_t i : range(nbin)) {
            auto last = std::partition(first, tmp.data.end(), [&bins,i](rtype t) {
                return t >= bins.safe(0,i) && t < bins.safe(1,i);
            });

            counts.safe[i] = last - first;
            first = last;

            if (last == tmp.data.end()) break;
        }

        return counts;
    }

    template<std::size_t Dim, typename Type, typename TypeB, typename TypeW>
    vec<1,meta::rtype_t<TypeW>> histogram(const vec<Dim,Type>& data, const vec<Dim,TypeW>& weight,
        const vec<2,TypeB>& bins) {
        phypp_check(bins.dims[0] == 2, "can only be called with a bin vector (expected "
            "dims=[2, ...], got dims=[", bins.dims, "])");
        phypp_check(data.dims == weight.dims, "incompatible dimensions for data and weight "
            "(", data.dims, " vs. ", weight.dims, ")");

        vec1u ids = uindgen(data.size());

        uint_t nbin = bins.dims[1];
        vec<1,meta::rtype_t<TypeW>> counts(nbin);

        auto first = ids.data.begin();
        for (uint_t i : range(nbin)) {
            auto last = std::partition(first, ids.data.end(), [&bins,&data,i](uint_t id) {
                return data.safe[id] >= bins.safe(0,i) && data.safe[id] < bins.safe(1,i);
            });

            for (; first != last; ++first) {
                counts.safe[i] += weight.safe[*first];
            }

            if (last == ids.data.end()) break;
        }

        return counts;
    }

    namespace impl {
        template<std::size_t Dim, typename Type, typename TypeB, typename F>
        void histogram_impl(const vec<Dim,Type>& data, const vec<2,TypeB>& bins, F&& func) {
            phypp_check(bins.dims[0] == 2, "can only be called with a bin vector (expected "
                "dims=[2, ...], got dims=[", bins.dims, "])");

            vec1u ids = uindgen(data.size());
            using iterator = vec1u::const_iterator;

            uint_t nbin = bins.dims[1];
            auto first = ids.data.begin();
            for (uint_t i : range(nbin)) {
                auto last = std::partition(first, ids.data.end(), [&bins,&data,i](uint_t id) {
                    return data.safe[id] >= bins.safe(0,i) && data.safe[id] < bins.safe(1,i);
                });

                func(i, meta::add_const(ids), iterator{first}, iterator{last});

                if (last == ids.data.end()) break;
            }
        }
    }

    template<std::size_t Dim, typename Type, typename TypeB, typename TypeF>
    void histogram(const vec<Dim,Type>& x, const vec<2,TypeB>& bins, TypeF&& func) {
        using iterator = vec1u::const_iterator;
        vec1u tids;

        impl::histogram_impl(x, bins,
            [&](uint_t i, const vec1u& ids, iterator i0, iterator i1) {
                uint_t npts = i1 - i0;
                uint_t k0 = i0 - ids.data.begin();
                tids.resize(npts);
                for (uint_t k : range(npts)) {
                    tids.safe[k] = ids.safe[k0 + k];
                }

                func(i, std::move(tids));
            }
        );
    }

    namespace impl {
        template<std::size_t Dim, typename TypeX, typename TypeY, typename TypeBX,
            typename TypeBY, typename TypeF>
        void histogram2d_impl(const vec<Dim,TypeX>& x, const vec<Dim,TypeY>& y,
            const vec<2,TypeBX>& xbins, const vec<2,TypeBY>& ybins, TypeF&& func) {

            phypp_check(xbins.dims[0] == 2, "can only be called with a bin vector (expected "
                "dims=[2, ...], got dims=[", xbins.dims, "])");
            phypp_check(ybins.dims[0] == 2, "can only be called with a bin vector (expected "
                "dims=[2, ...], got dims=[", ybins.dims, "])");
            phypp_check(x.dims == y.dims, "incompatible dimensions for x and y (", x.dims, " vs. ",
                y.dims, ")");

            vec1u ids = uindgen(x.size());

            uint_t nxbin = xbins.dims[1];
            uint_t nybin = ybins.dims[1];

            using iterator = vec1u::const_iterator;

            auto firstx = ids.data.begin();
            for (uint_t i : range(nxbin)) {
                auto lastx = std::partition(firstx, ids.data.end(), [&xbins,&x,i](uint_t id) {
                    return x.safe[id] >= xbins.safe(0,i) && x.safe[id] < xbins.safe(1,i);
                });

                auto firsty = firstx;
                for (uint_t j : range(nybin)) {
                    auto lasty = std::partition(firsty, lastx, [&ybins,&y,j](uint_t id) {
                        return y.safe[id] >= ybins.safe(0,j) && y.safe[id] < ybins.safe(1,j);
                    });

                    func(i, j, meta::add_const(ids), iterator{firsty}, iterator{lasty});

                    if (lasty == lastx) break;
                    firsty = lasty;
                }

                if (lastx == ids.data.end()) break;
                firstx = lastx;
            }
        }
    }

    template<std::size_t Dim, typename TypeX, typename TypeY, typename TypeBX, typename TypeBY>
    vec2u histogram2d(const vec<Dim,TypeX>& x, const vec<Dim,TypeY>& y,
        const vec<2,TypeBX>& xbins, const vec<2,TypeBY>& ybins) {

        vec2u counts(xbins.dims[1], ybins.dims[1]);

        using iterator = vec1u::const_iterator;
        impl::histogram2d_impl(x, y, xbins, ybins,
            [&](uint_t i, uint_t j, const vec1u& ids, iterator i0, iterator i1) {
                counts.safe(i,j) = i1 - i0;
            }
        );

        return counts;
    }

    template<std::size_t Dim, typename TypeX, typename TypeY, typename TypeBX,
        typename TypeBY, typename TypeF>
    void histogram2d(const vec<Dim,TypeX>& x, const vec<Dim,TypeY>& y,
        const vec<2,TypeBX>& xbins, const vec<2,TypeBY>& ybins, TypeF&& func) {

        using iterator = vec1u::const_iterator;
        impl::histogram2d_impl(x, y, xbins, ybins,
            [&](uint_t i, uint_t j, const vec1u& ids, iterator i0, iterator i1) {
                uint_t npts = i1 - i0;
                vec1u tids(npts);
                uint_t k0 = i0 - ids.data.begin();
                for (uint_t k : range(npts)) {
                    tids.safe[k] = ids.safe[k0 + k];
                }

                func(i, j, std::move(tids));
            }
        );
    }

    template<typename TypeY>
    auto cumul(const vec<1,TypeY>& y) -> vec<1,meta::rtype_t<TypeY>> {
        vec<1,meta::rtype_t<TypeY>> dr(y.dims);
        for (uint_t i : range(1, y.size())) {
            dr.safe[i] = dr.safe[i-1] + y.safe[i];
        }

        return dr;
    }

    template<typename TypeY>
    auto cumul_reverse(const vec<1,TypeY>& y) -> vec<1,meta::rtype_t<TypeY>> {
        vec<1,meta::rtype_t<TypeY>> dr(y.dims);
        for (uint_t i : range(1, y.size())) {
            dr.safe[dr.size()-1-i] = dr.safe[dr.size()-i] + y.safe[y.size()-1-i];
        }

        return dr;
    }
}

#endif
