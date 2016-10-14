#ifndef PHYPP_MATH_REDUCE_HPP
#define PHYPP_MATH_REDUCE_HPP

#include "phypp/core/vec.hpp"
#include "phypp/core/range.hpp"
#include "phypp/core/error.hpp"
#include "phypp/core/string_conversion.hpp"
#include "phypp/utility/generic.hpp"
#include "phypp/math/base.hpp"

namespace phypp {
    namespace meta {
        template<typename T>
        using total_return_type = typename std::conditional<std::is_integral<T>::value,
            typename std::conditional<std::is_unsigned<T>::value, uint_t, int_t>::type,
            double>::type;
    }

    template<std::size_t Dim, typename Type>
    meta::total_return_type<meta::rtype_t<Type>> total(const vec<Dim,Type>& v) {
        meta::total_return_type<meta::rtype_t<Type>> total = 0;
        for (auto& t : v) {
            total += t;
        }

        return total;
    }

    template<std::size_t Dim = 1, typename Type = bool, typename enable =
        typename std::enable_if<std::is_same<meta::rtype_t<Type>, bool>::value>::type>
    uint_t count(const vec<Dim,Type>& v) {
        uint_t n = 0u;
        for (bool b : v) {
            if (b) ++n;
        }

        return n;
    }

    template<std::size_t Dim, typename Type>
    double mean(const vec<Dim,Type>& v) {
        double total = 0.0;
        for (auto& t : v) {
            total += t;
        }

        return total/v.size();
    }

    template<std::size_t Dim, typename Type, typename TypeW>
    double weighted_mean(const vec<Dim,Type>& v, const vec<Dim,TypeW>& w) {
        phypp_check(v.dims == w.dims, "incompatible dimensions between values and weights "
            "(", v.dims, " vs. ", w.dims, ")");

        double total = 0.0;
        double totalw = 0.0;
        for (uint_t i : range(v)) {
            total += v.safe[i]*w.safe[i];
            totalw += w.safe[i];
        }

        return total/totalw;
    }

    template<std::size_t Dim, typename Type, typename TypeW>
    std::pair<double,double> optimal_mean(const vec<Dim,Type>& v, const vec<Dim,TypeW>& e) {
        phypp_check(v.dims == e.dims, "incompatible dimensions between values and uncertianties "
            "(", v.dims, " vs. ", e.dims, ")");

        double totalv = 0.0;
        double totalw = 0.0;
        for (uint_t i : range(v)) {
            double w = 1.0/(e.safe[i]*e.safe[i]);
            totalv += w*v.safe[i];
            totalw += w;
        }

        return std::make_pair(totalv/totalw, 1/sqrt(totalw));
    }

    template<std::size_t Dim, typename T, typename enable =
        typename std::enable_if<std::is_same<meta::rtype_t<T>, bool>::value>::type>
    double fraction_of(const vec<Dim,T>& b) {
        return mean(b);
    }

    template<std::size_t Dim, typename T, typename U, typename enable =
        typename std::enable_if<std::is_same<meta::rtype_t<T>, bool>::value &&
                                std::is_same<meta::rtype_t<U>, bool>::value>::type>
    double fraction_of(const vec<Dim,T>& b, const vec<Dim,U>& among) {
        phypp_check(b.dims == among.dims, "incompatible dimensions between count vector and "
            "among vector (", b.dims, " vs. ", among.dims, ")");

        double total = 0.0;
        uint_t count = 0;
        for (uint_t i : range(among)) {
            if (among.safe[i]) {
                total += b.safe[i];
                ++count;
            }
        }

        return total/count;
    }

    template<typename ... Args>
    std::string percent_of(Args&& ... args) {
        return strn(round(100.0*fraction_of(std::forward<Args>(args)...)))+"%";
    }

    template<std::size_t Dim, typename Type>
    meta::rtype_t<Type> inplace_median(vec<Dim,Type>& v) {
        phypp_check(!v.empty(), "cannot find the median of an empty vector");

        using vtype = typename vec<Dim,Type>::vtype;
        using dtype = typename vtype::value_type;

        uint_t nwrong = 0;
        for (auto& t : v) {
            nwrong += is_nan(t);
        }

        if (nwrong == v.size()) return dnan;

        std::ptrdiff_t offset = (v.size()-nwrong)/2;
        std::nth_element(v.data.begin(), v.data.begin() + offset, v.data.end(),
            [](dtype i, dtype j) {
                if (is_nan(impl::dref<Type>(i))) return false;
                if (is_nan(impl::dref<Type>(j))) return true;
                return impl::dref<Type>(i) < impl::dref<Type>(j);
            }
        );

        return *(v.begin() + offset);
    }

    template<std::size_t Dim, typename Type>
    meta::rtype_t<Type> median(vec<Dim,Type> v) {
        return inplace_median(v);
    }

    template<std::size_t Dim, typename Type, typename TypeW>
    meta::rtype_t<Type> weighted_median(const vec<Dim,Type>& v, const vec<Dim,TypeW>& w) {
        phypp_check(!v.empty(), "cannot find the weighted median of an empty vector");
        phypp_check(v.dims == w.dims, "incompatible dimensions between values and weights "
            "(", v.dims, " vs. ", w.dims, ")");

        double totw = 0;
        for (uint_t i : range(v)) {
            if (!is_nan(v.safe[i]) && !is_nan(w.safe[i])) totw += w.safe[i];
        }

        vec1u ids = sort(v);

        double tot = 0;
        for (uint_t i : range(ids)) {
            uint_t j = ids.safe[i];
            if (!is_nan(v.safe[j]) && !is_nan(w.safe[j])) {
                tot += w.safe[j];
                if (tot > totw/2.0) {
                    if (i == 0) return v.safe[j];
                    else        return v.safe[ids.safe[i-1]];
                }
            }
        }

        return dnan;
    }

    template<std::size_t Dim, typename Type, typename U>
    meta::rtype_t<Type> percentile(const vec<Dim,Type>& v, const U& u) {
        phypp_check(!v.empty(), "cannot find the percentiles of an empty vector");

        vec1u ok = where(is_finite(v));
        if (ok.empty()) return 0;

        // TODO: (fixme) use same algorithm than median
        typename vec<1,Type>::effective_type t = v.safe[ok];
        std::ptrdiff_t offset = clamp(t.size()*u, 0u, t.size()-1);
        std::nth_element(t.begin(), t.begin() + offset, t.end());
        return *(t.begin() + offset);
    }

    namespace impl {
        template<std::size_t Dim, typename Type>
        void percentiles_(vec<1,Type>& r, uint_t i, vec<Dim,Type>& t) {}

        template<std::size_t Dim, typename Type, typename U, typename ... Args>
        void percentiles_(vec<1,Type>& r, uint_t i, vec<Dim,Type>& t, const U& u, const Args& ... args) {
            std::ptrdiff_t offset = clamp(t.size()*u, 0u, t.size()-1);
            std::nth_element(t.begin(), t.begin() + offset, t.end());
            r.safe[i] = *(t.begin() + offset);
            ++i;

            percentiles_(r, i, t, args...);
        }
    }

    template<std::size_t Dim, typename Type, typename ... Args>
    typename vec<1,Type>::effective_type percentiles(const vec<Dim,Type>& v, const Args& ... args) {
        phypp_check(!v.empty(), "cannot find the percentiles of an empty vector");

        vec1u ok = where(is_finite(v));
        typename vec<1,Type>::effective_type t;
        if (ok.empty()) return t;
        t = v.safe[ok];

        typename vec<1,Type>::effective_type r = arr<meta::rtype_t<Type>>(sizeof...(Args));
        impl::percentiles_(r, 0, t, args...);

        return r;
    }

    template<std::size_t Dim, typename Type>
    vec<Dim,bool> sigma_clip(const vec<Dim,Type>& tv, double sigma) {
        auto v = tv.concretise();
        auto med = inplace_median(v);
        auto mad = median(abs(v - med));
        // Note: cannot use 'v' below, since the order of the values has changed!
        return abs(tv - med) < sigma*mad;
    }

    namespace impl {
        template<std::size_t Dim, typename Type>
        typename vec<Dim,Type>::const_iterator min_(const vec<Dim,Type>& v) {
            phypp_check(!v.empty(), "cannot find the minimum of an empty vector");

            auto iter = std::min_element(v.begin(), v.end(), [](meta::rtype_t<Type> t1, meta::rtype_t<Type> t2){
                if (is_nan(t1)) return false;
                if (is_nan(t2)) return true;
                return t1 < t2;
            });

            if (iter == v.end()) iter = v.begin();
            return iter;
        }

        template<std::size_t Dim, typename Type>
        typename vec<Dim,Type>::const_iterator max_(const vec<Dim,Type>& v) {
            phypp_check(!v.empty(), "cannot find the maximum of an empty vector");

            auto iter = std::max_element(v.begin(), v.end(), [](meta::rtype_t<Type> t1, meta::rtype_t<Type> t2){
                if (is_nan(t1)) return true;
                if (is_nan(t2)) return false;
                return t1 < t2;
            });

            if (iter == v.end()) iter = v.begin();
            return iter;
        }

        template<std::size_t Dim, typename Type>
        std::pair<typename vec<Dim,Type>::const_iterator, typename vec<Dim,Type>::const_iterator>
            minmax_(const vec<Dim,Type>& v) {
            phypp_check(!v.empty(), "cannot find the maximum/minimum of an empty vector");

            // We cannot take care of NaN using std::minmax_element and the trick of
            // std::min_element and std::max_element. So we just roll our own...
            // This naive version performs slightly more comparisons than the standard
            // algorithm. This should not be noticeable, but could be improved by using a
            // filtered iterator that automatically skips NaN values.
            // See, e.g., boost::filter_iterator for a working implementation.
            using iterator = typename vec<Dim,Type>::const_iterator;
            std::pair<iterator,iterator> res;

            iterator i = v.begin();
            while (i != v.end() && is_nan(*i)) {
                ++i;
            }

            if (i == v.end()) {
                // Only NaN
                res.first = v.begin(); res.second = v.begin();
                return res;
            }

            res.first = i; res.second = i;

            for (++i; i != v.end(); ++i) {
                if (is_nan(*i)) continue;
                if (*i < *res.first)          res.first = i;
                else if (!(*i < *res.second)) res.second = i;
            }

            return res;
        }
    }

    template<std::size_t Dim, typename Type>
    meta::rtype_t<Type> min(const vec<Dim,Type>& v) {
        return *impl::min_(v);
    }

    template<std::size_t Dim, typename Type>
    meta::rtype_t<Type> max(const vec<Dim,Type>& v) {
        return *impl::max_(v);
    }

    template<std::size_t Dim, typename Type>
    std::pair<meta::rtype_t<Type>,meta::rtype_t<Type>> minmax(const vec<Dim,Type>& v) {
        auto tmp = impl::minmax_(v);
        return std::make_pair(*tmp.first, *tmp.second);
    }

    template<std::size_t Dim, typename Type>
    meta::rtype_t<Type> min(const vec<Dim,Type>& v, uint_t& id) {
        auto iter = impl::min_(v);
        id = iter - v.begin();
        return *iter;
    }

    template<std::size_t Dim, typename Type>
    meta::rtype_t<Type> max(const vec<Dim,Type>& v, uint_t& id) {
        auto iter = impl::max_(v);
        id = iter - v.begin();
        return *iter;
    }

    template<std::size_t Dim, typename Type>
    std::pair<meta::rtype_t<Type>,meta::rtype_t<Type>> minmax(const vec<Dim,Type>& v, std::pair<uint_t,uint_t>& ids) {
        auto tmp = impl::minmax_(v);
        ids.first = tmp.first - v.begin();
        ids.second = tmp.second - v.begin();
        return std::make_pair(*tmp.first, *tmp.second);
    }

    template<std::size_t Dim, typename Type>
    uint_t min_id(const vec<Dim,Type>& v) {
        return impl::min_(v) - v.begin();
    }

    template<std::size_t Dim, typename Type>
    uint_t max_id(const vec<Dim,Type>& v) {
        return impl::max_(v) - v.begin();
    }

    template<std::size_t Dim, typename Type>
    std::pair<uint_t,uint_t> minmax_ids(const vec<Dim,Type>& v) {
        auto tmp = impl::minmax_(v);
        return std::make_pair(uint_t(tmp.first - v.begin()), uint_t(tmp.second - v.begin()));
    }

    template<typename T1, typename T2>
    auto min(T1&& t1, T2&& t2) -> decltype((t1 < t2 ? t1 : t2)) {
        return t1 < t2 ? t1 : t2;
    }

    template<typename T1, typename T2>
    auto max(T1&& t1, T2&& t2) -> decltype((t1 > t2 ? t1 : t2)) {
        return t1 > t2 ? t1 : t2;
    }

    template<std::size_t Dim, typename Type1, typename Type2>
    vec<Dim,meta::rtype_t<Type1>> min(const vec<Dim,Type1>& v1, const vec<Dim,Type2>& v2) {
        phypp_check(v1.dims == v2.dims, "min: incompatible vector dimensions "
            "(", v1.dims, " vs. ", v2.dims, ")");

        vec<Dim,meta::rtype_t<Type1>> r(v1.dims);
        for (uint_t i = 0; i < v1.size(); ++i) {
            r.safe[i] = std::min<decltype(v1[0]*v2[0])>(v1.safe[i], v2.safe[i]);
        }
        return r;
    }

    template<std::size_t Dim, typename Type1, typename Type2>
    vec<Dim,meta::rtype_t<Type1>> max(const vec<Dim,Type1>& v1, const vec<Dim,Type2>& v2) {
        phypp_check(v1.dims == v2.dims, "max: incompatible vector dimensions "
            "(", v1.dims, " vs. ", v2.dims, ")");

        vec<Dim,meta::rtype_t<Type1>> r(v1.dims);
        for (uint_t i = 0; i < v1.size(); ++i) {
            r.safe[i] = std::max<decltype(v1[0]*v2[0])>(v1.safe[i], v2.safe[i]);
        }
        return r;
    }

    template<std::size_t Dim, typename Type1, typename Type2>
    vec<Dim,meta::rtype_t<Type1>> min(const vec<Dim,Type1>& v1, const Type2& v2) {
        vec<Dim,meta::rtype_t<Type1>> r(v1.dims);
        for (uint_t i = 0; i < v1.size(); ++i) {
            r.safe[i] = std::min<decltype(v1[0]*v2)>(v1.safe[i], v2);
        }
        return r;
    }

    template<std::size_t Dim, typename Type1, typename Type2>
    vec<Dim,meta::rtype_t<Type1>> max(const vec<Dim,Type1>& v1, const Type2& v2) {
        vec<Dim,meta::rtype_t<Type1>> r(v1.dims);
        for (uint_t i = 0; i < v1.size(); ++i) {
            r.safe[i] = std::max<decltype(v1[0]*v2)>(v1.safe[i], v2);
        }
        return r;
    }

    template<std::size_t Dim, typename Type1, typename Type2>
    vec<Dim,meta::rtype_t<Type1>> min(const Type1& v1, const vec<Dim,Type2>& v2) {
        vec<Dim,meta::rtype_t<Type1>> r(v2.dims);
        for (uint_t i = 0; i < v2.size(); ++i) {
            r.safe[i] = std::min<decltype(v1*v2[0])>(v1, v2.safe[i]);
        }
        return r;
    }

    template<std::size_t Dim, typename Type1, typename Type2>
    vec<Dim,meta::rtype_t<Type1>> max(const Type1& v1, const vec<Dim,Type2>& v2) {
        vec<Dim,meta::rtype_t<Type1>> r(v2.dims);
        for (uint_t i = 0; i < v2.size(); ++i) {
            r.safe[i] = std::max<decltype(v1*v2[0])>(v1, v2.safe[i]);
        }
        return r;
    }

    template<std::size_t Dim, typename Type>
    double rms(const vec<Dim,Type>& v) {
        double sum = 0;
        for (auto& t : v) {
            sum += t*t;
        }

        return sqrt(sum/v.size());
    }

    template<std::size_t Dim, typename Type>
    double stddev(const vec<Dim,Type>& v) {
        return rms(v - mean(v));
    }

    template<std::size_t Dim, typename Type>
    meta::rtype_t<Type> mad(const vec<Dim,Type>& v) {
        return median(abs(v - median(v)));
    }

    namespace impl {
        template<typename F, F f, std::size_t Dim, typename Type, typename ... Args>
        auto run_index_(uint_t dim, const vec<Dim,Type>& v, Args&& ... args) ->
        vec<Dim-1,typename meta::return_type<F>::type> {

            vec<Dim-1,typename meta::return_type<F>::type> r;
            for (uint_t i = 0; i < dim; ++i) {
                r.dims[i] = v.dims[i];
            }
            for (uint_t i = dim+1; i < Dim; ++i) {
                r.dims[i-1] = v.dims[i];
            }
            r.resize();

            uint_t mpitch = 1;
            for (uint_t i = dim+1; i < Dim; ++i) {
                mpitch *= v.dims[i];
            }

        //  Example demonstration of the index computation:
        //
        //  Assume we have a 4D array of dimensions d1, d2, d3 and d4.
        //  We want to compute an operation on the array generated by:
        //      data[i*d2*d3*d4 + j*d3*d4 + k*d4 + l]
        //  and by running the 'dim'th index (i, j, k or l), assigning
        //  the value to another 3D array using a continuous index 'u'.
        //
        //  - (i,j,k,_) : k+j*d3+i*d2*d3 = u : mean(((u/1)*d4 + l)*1 + (u%1))
        //      mpitch = 1
        //
        //  - (i,j,_,l) : l+j*d4+i*d2*d4 = u : mean(((u/d4)*d3 + k)*d4 + (u%d4))
        //      mpitch = d4
        //
        //  - (i,_,k,l) : l+k*d4+i*d3*d4 = u : mean(((u/(d3*d4))*d2 + j)*d3*d4 + (u%(d3*d4)))
        //      mpitch = d4*d3
        //
        //  - (_,j,k,l) : l+k*d4+j*d3*d4 = u : mean(((u/(d4*d3*d2))*d1 + i)*d2*d3*d4 + (u%(d4*d3*d2)))
        //      mpitch = d4*d3*d2
        //
        //  Final recipe:
        //      ((u/mpitch)*dim[d] + i)*mpitch + (u%mpitch)

            vec<1,meta::rtype_t<Type>> tmp(v.dims[dim]);
            for (uint_t i : range(r)) {
                uint_t base = (i%mpitch) + (i/mpitch)*v.dims[dim]*mpitch;
                for (uint_t j : range(tmp)) {
                    tmp.safe[j] = v.safe[base + j*mpitch];
                }

                r.safe[i] = (*f)(tmp, std::forward<Args>(args)...);
            }

            return r;
        }

        template<std::size_t Dim, typename Type>
        vec<1,meta::rtype_t<Type>> run_dim_apply_ids_(const vec1u& ids, const vec<Dim,Type>& v) {
            return v.safe[ids].concretise();
        }

        template<std::size_t Dim, typename Type, typename ... Args>
        std::array<uint_t,Dim> run_dim_get_dim_(const vec<Dim,Type>& v, const Args& ... vs) {
            return v.dims;
        }

        template<typename F, typename ... Args>
        void run_dim_final_(uint_t dim, F&& func, const Args& ... vs) {
            auto ds = run_dim_get_dim_(vs...);
            const uint_t N = meta::array_size<decltype(ds)>::size;

            uint_t nint = ds[dim];
            uint_t mpitch = 1;
            uint_t np = 1;
            for (uint_t i : range(N)) {
                if (i != dim) np *= ds[i];
                if (i > dim) mpitch *= ds[i];
            }

            vec1u ids(nint);
            for (uint_t i : range(np)) {
                uint_t base = (i%mpitch) + (i/mpitch)*nint*mpitch;
                for (uint_t j : range(ids)) {
                    ids.safe[j] = base + j*mpitch;
                }

                func(i, run_dim_apply_ids_<N>(ids, vs)...);
            }
        }

        template<typename ... Args1>
        struct run_dim_unroll_ {
            template<typename T>
            static void run(uint_t dim, Args1&&... a1, impl::placeholder_t, T&& t) {
                run_dim_final_(dim, std::forward<T>(t), std::forward<Args1>(a1)...);
            }

            template<typename T, typename ... Args2>
            static void run(uint_t dim, Args1&&... a1, impl::placeholder_t, T&& t, Args2&&... a2) {
                run_dim_unroll_<Args1..., T>::run(dim, std::forward<Args1>(a1)...,
                    std::forward<T>(t), _, std::forward<Args2>(a2)...);
            }
        };

        template<std::size_t D1, typename F>
        bool run_dim_check_dims__(const std::array<uint_t,D1>& dims, F&& func) {
            return true;
        }

        template<std::size_t D1, std::size_t D2, typename Type>
        bool run_dim_check_dims__(const std::array<uint_t,D1>& dims, const vec<D2,Type>& v) {
            return dims == v.dims;
        }

        template<std::size_t Dim, typename Type, typename ... Args>
        bool run_dim_check_dims_(uint_t& dim, const vec<Dim,Type>& v, const Args& ... vs) {
            dim = Dim;
            return count({!run_dim_check_dims__(v.dims, vs)...}) == 0;
        }
    }

    // Iterate over one dimension of the provided vector and call a function for each slice.
    template<typename ... Args>
    void run_dim(uint_t dim, Args&& ... args) {
        uint_t Dim;
        bool check = impl::run_dim_check_dims_(Dim, args...);
        phypp_check(check, "incompatible dimensions of input vectors");
        phypp_check(dim < Dim, "reduction dimension is incompatible with input vectors");

        impl::run_dim_unroll_<>::run(dim, _, std::forward<Args>(args)...);
    }

    template<typename F, std::size_t Dim, typename Type>
    auto reduce(uint_t dim, const vec<Dim,Type>& v, F&& func) ->
        vec<Dim-1,typename meta::return_type<F>::type> {
        phypp_check(dim < Dim, "reduction dimension is incompatible with input vector "
            "(", dim, " vs. ", v.dims, ")");

        vec<Dim-1,typename meta::return_type<F>::type> r;
        for (uint_t i = 0; i < dim; ++i) {
            r.dims[i] = v.dims[i];
        }
        for (uint_t i = dim+1; i < Dim; ++i) {
            r.dims[i-1] = v.dims[i];
        }

        r.resize();

        run_dim(dim, v, [&](uint_t i, vec<1,meta::rtype_t<Type>> tv) {
            r.safe[i] = func(std::move(tv));
        });

        return r;
    }

    #define MAKE_PARTIAL(func) \
        template<typename T, typename ... Args> \
        struct func ## _run_index_wrapper_ { \
            static auto run(const vec<1,T>& v, Args&& ... args) -> \
            decltype(func(v, std::forward<Args>(args)...)) { \
                return func(v, std::forward<Args>(args)...); \
            } \
        }; \
        \
        template<std::size_t Dim, typename Type, typename ... Args> \
        auto partial_ ## func (uint_t dim, const vec<Dim,Type>& v, Args&& ... args) -> \
        vec<Dim-1, decltype(func(std::declval<vec<1,meta::rtype_t<Type>>>(), std::forward<Args>(args)...))> { \
            phypp_check(dim < Dim, "reduction dimension is incompatible with input vector " \
                "(", dim, " vs. ", v.dims, ")"); \
            using wrapper = func ## _run_index_wrapper_<meta::rtype_t<Type>, Args...>; \
            using fptr = decltype(&wrapper::run); \
            return impl::run_index_<fptr, &wrapper::run>(dim, v, std::forward<Args>(args)...); \
        }

    MAKE_PARTIAL(total);
    MAKE_PARTIAL(count);
    MAKE_PARTIAL(fraction_of);
    MAKE_PARTIAL(mean);
    MAKE_PARTIAL(median);
    MAKE_PARTIAL(min);
    MAKE_PARTIAL(max);
    MAKE_PARTIAL(percentile);
    MAKE_PARTIAL(rms);
    MAKE_PARTIAL(stddev);
    MAKE_PARTIAL(mad);

    #undef MAKE_PARTIAL
    
    namespace impl {
        template<std::size_t Dim, typename Type>
        void data_info_(const vec<Dim,Type>& v) {
            vec1u idok = where(is_finite(v));
            print(idok.size(), "/", v.size(), " valid values (dims: ", v.dims, ")");
            if (idok.size() == 0) return;

            vec<1,meta::rtype_t<Type>> tv = v.safe[idok];

            print(" min : ", min(tv));
            print(" 15% : ", percentile(tv, 0.15));
            print(" 50% : ", median(tv));
            print(" mean: ", mean(tv));
            print(" 85% : ", percentile(tv, 0.85));
            print(" max : ", max(tv));
            print(" rms : ", stddev(tv - median(tv)));
        }
    }

    #define data_info(x) \
        print("data info: ", #x, " (", typeid(x).name(), ")"); \
        phypp::impl::data_info_(x);

    template<typename TypeX, typename TypeY>
    auto integrate(const vec<1,TypeX>& x, const vec<1,TypeY>& y)
        -> decltype(0.5*y[0]*(x[1]-x[0])) {

        phypp_check(x.size() == y.size(),
            "incompatible x and y array dimensions (", x.size(), " vs ", y.size(), ")");

        decltype(0.5*y[0]*(x[1]-x[0])) r = 0;
        for (uint_t i = 0; i < x.size()-1; ++i) {
            r += 0.5*(y.safe[i+1]+y.safe[i])*(x.safe[i+1]-x.safe[i]);
        }

        return r;
    }

    template<typename TypeX, typename TypeY>
    auto integrate(const vec<2,TypeX>& x, const vec<1,TypeY>& y)
        -> decltype(y[0]*(x[1]-x[0])) {

        phypp_check(x.dims[0] == 2, "x array must be a binned array (expected dim=[2,...], "
            "got dim=[", x.dims[0], ",...]");
        phypp_check(x.dims[1] == y.size(),
            "incompatible x and y array dimensions (", x.dims[1], " vs ", y.size(), ")");

        decltype(y[0]*(x[1]-x[0])) r = 0;
        for (uint_t i = 0; i < x.dims[1]; ++i) {
            r += y.safe[i]*(x.safe[i+x.dims[1]]-x.safe[i]);
        }

        return r;
    }

    template<typename TypeX, typename TypeY>
    auto integrate(const vec<1,TypeX>& x, const vec<1,TypeY>& y, double x0, double x1)
        -> decltype(0.5*y[0]*(x[1]-x[0])) {

        phypp_check(x.size() == y.size(),
            "incompatible x and y array dimensions (", x.size(), " vs ", y.size(), ")");

        phypp_check(x.front() <= x0 && x.back() >= x1, "x array must cover the range [x0,x1]");

        uint_t i0 = upper_bound(x0, x);
        uint_t i1 = lower_bound(x1, x);

        if (i0 > i1) {
            return 0.5*(y.safe[i1]+y.safe[i0])*(x1 - x0);
        } else {
            decltype(0.5*y[0]*(x[1]-x[0])) r = 0;
            for (uint_t i = i0; i < i1; ++i) {
                r += 0.5*(y.safe[i+1]+y.safe[i])*(x.safe[i+1]-x.safe[i]);
            }

            if (i0 > 0) {
                double y0 = interpolate(y.safe[i0-1], y.safe[i0], x.safe[i0-1], x.safe[i0], x0);
                r += 0.5*(y.safe[i0]+y0)*(x[i0]-x0);
            }

            if (i1 < x.size()-1) {
                double y1 = interpolate(y.safe[i1], y.safe[i1+1], x.safe[i1], x.safe[i1+1], x1);
                r += 0.5*(y1+y.safe[i1])*(x1-x.safe[i1]);
            }

            return r;
        }
    }

    template<typename TypeX, typename TypeY>
    auto integrate(const vec<2,TypeX>& x, const vec<1,TypeY>& y, double x0, double x1)
        -> decltype(y[0]*(x[1]-x[0])) {

        phypp_check(x.dims[0] == 2, "x array must be a binned array (expected dim=[2,...], "
            "got dim=[", x.dims[0], ",...]");
        phypp_check(x.dims[1] == y.size(),
            "incompatible x and y array dimensions (", x.dims[1], " vs ", y.size(), ")");

        phypp_check(x.safe[0] <= x0 && x.safe[x.size()-1] >= x1, "x array must cover the range [x0,x1]");

        uint_t i0 = upper_bound(x0, x(0,_));
        uint_t i1 = lower_bound(x1, x(0,_));

        if (i0 > i1) {
            return y.safe[i1]*(x1 - x0);
        } else {
            decltype(y[0]*(x[1]-x[0])) r = 0;
            for (uint_t i = i0; i < i1; ++i) {
                r += y.safe[i]*(x.safe[i+x.dims[1]]-x.safe[i]);
            }

            r += y.safe[i0]*(x.safe[i0] - x0);
            r += y.safe[i1]*(x1 - x.safe[i1]);

            return r;
        }
    }

    template<typename TypeX, typename TypeY>
    auto cumul(const vec<1,TypeX>& x, const vec<1,TypeY>& y) -> vec<1,decltype(integrate(x,y))> {
        vec<1,decltype(integrate(x,y))> dr(y.dims);

        phypp_check(x.size() == y.size(),
            "incompatible x and y array dimensions (", x.size(), " vs ", y.size(), ")");

        decltype(0.5*y[0]*(x[1]-x[0])) r = 0;
        for (uint_t i : range(1, x.size())) {
            r += 0.5*(y.safe[i]+y.safe[i-1])*(x.safe[i]-x.safe[i-1]);
            dr.safe[i] = r;
        }

        return dr;
    }

    template<typename F, typename T, typename U,
        typename enable = typename std::enable_if<!meta::is_vec<F>::value>::type>
    auto integrate_func(F f, T x0, U x1,
        double e = std::numeric_limits<typename std::result_of<F(T)>::type>::epsilon()) ->
        decltype(0.5*f(x0)*x0) {

        using rtype = decltype(0.5*f(x0)*x0);

        std::vector<rtype> buffer;
        buffer.reserve(20);
        buffer.push_back(0.5*(x1 - x0)*(f(x0) + f(x1)));

        uint_t n = 0;
        uint_t oid = 0;

        do {
            ++n;

            rtype tr = 0;
            uint_t tn = 1 << n;
            auto d = (x1 - x0)/double(tn);
            for (uint_t k = 1; k <= tn/2; ++k) {
                tr += f(x0 + (2*k-1)*d);
            }
            buffer.push_back(0.5*buffer[oid] + d*tr);

            for (uint_t m = 1; m <= n; ++m) {
                auto t = 1 << (2*m);
                buffer.push_back((t*buffer.back() - buffer[oid+m-1])/(t - 1));
            }

            oid += n;
        } while (abs((buffer.back() - *(buffer.end()-2))/buffer.back()) > e);

        return buffer.back();
    }
}

#endif
