#ifndef PHYPP_INCLUDING_GENERIC_BITS
#error this file is not meant to be included separately, include "phypp/utilty/generic.hpp" instead
#endif

namespace phypp {
    template<typename Type>
    vec<1,Type> reverse(vec<1,Type> v) {
        std::reverse(v.data.begin(), v.data.end());
        return v;
    }

    template<typename Type,
        typename enable = typename std::enable_if<!std::is_pointer<Type>::value>::type>
    void inplace_shift(vec<1,Type>& v, int_t n) {
        n = (-n) % int_t(v.size());
        if (n < 0) n = int_t(v.size())+n;

        if (n != 0) {
            std::rotate(v.data.begin(), v.data.begin() + n, v.end());
        }
    }

    template<typename Type,
        typename enable = typename std::enable_if<!std::is_pointer<Type>::value>::type>
    vec<1,Type> shift(vec<1,Type> v, int_t n) {
        inplace_shift(v, n);
        return v;
    }

    template<typename Type,
        typename enable = typename std::enable_if<!std::is_pointer<Type>::value>::type>
    vec<2,Type> transpose(const vec<2,Type>& v) {
        vec<2,Type> r(impl::vec_nocopy_tag, v);
        std::swap(r.dims[0], r.dims[1]);

        for (uint_t i : range(v)) {
            r.data.push_back(v.data[(i%v.dims[0])*v.dims[1] + i/v.dims[0]]);
        }

        // TODO: (optimization) see who's faster
        // r.resize();
        // for (uint_t i : range(r.dims[0]))
        // for (uint_t j : range(r.dims[1])) {
        //     r.data[j+i*r.dims[1]] = v.data[i+j*v.dims[1]];
        // }

        return r;
    }

    template<std::size_t Dim, typename Type = double, typename ... Args>
    vec<Dim+meta::dim_total<Args...>::value, meta::rtype_t<Type>>
        replicate(const vec<Dim,Type>& t, Args&& ... args) {
        static const std::size_t FDim = Dim+meta::dim_total<Args...>::value;
        vec<FDim, meta::rtype_t<Type>> v(std::forward<Args>(args)..., t.dims);

        std::size_t pitch = t.size();
        std::size_t n = v.size()/pitch;
        for (uint_t i : range(n))
        for (uint_t j : range(pitch)) {
            v.safe[i*pitch + j] = t.safe[j];
        }

        return v;
    }

    template<typename Type, typename ... Args>
    vec<meta::dim_total<Args...>::value, meta::vtype_t<Type>> replicate(const Type& t, Args&& ... args) {
        static const std::size_t FDim = meta::dim_total<Args...>::value;
        vec<FDim, meta::vtype_t<Type>> v(std::forward<Args>(args)...);

        for (auto& e : v) {
            e = t;
        }

        return v;
    }

    template<std::size_t Dim, typename Type>
    vec1u sort(const vec<Dim,Type>& v) {
        vec1u r = uindgen(v.size());
        std::stable_sort(r.data.begin(), r.data.end(), [&v](uint_t i, uint_t j) {
            return typename vec<Dim,Type>::comparator_less()(v.data[i], v.data[j]);
        });

        return r;
    }

    template<std::size_t Dim, typename Type, typename F>
    vec1u sort(const vec<Dim,Type>& v, F&& comp) {
        vec1u r = uindgen(v.size());
        std::stable_sort(r.data.begin(), r.data.end(), [&v,&comp](uint_t i, uint_t j) {
            return comp(v.safe[i], v.safe[j]);
        });

        return r;
    }

    template<std::size_t Dim, typename Type>
    void inplace_sort(vec<Dim,Type>& v) {
        std::stable_sort(v.data.begin(), v.data.end(), typename vec<Dim,Type>::comparator_less());
    }

    template<std::size_t Dim, typename Type, typename F>
    void inplace_sort(vec<Dim,Type>& v, F&& comp) {
        std::stable_sort(v.begin(), v.end(), std::forward<F>(comp));
    }

    // Check if a given array is sorted or not
    template<std::size_t Dim, typename Type>
    bool is_sorted(const vec<Dim,Type>& v) {
        for (uint_t i : range(1, v.size())) {
            if (!typename vec<Dim,Type>::comparator_less()(v.safe[i-1], v.safe[i])) return false;
        }

        return true;
    }

    template<std::size_t Dim, typename Type,
        typename enable = typename std::enable_if<!std::is_pointer<Type>::value>::type>
    void inplace_remove(vec<Dim,Type>& v, vec1u ids) {
        inplace_sort(ids);

        uint_t i = 0;
        uint_t pitch = v.pitch(0);
        uint_t osize = v.size();
        while (i < ids.size()) {
            uint_t i1 = ids.safe[ids.size()-1-i];
            uint_t i0 = i1;

            phypp_check(i1 < osize, "trying to erase value ", i1, " in vector of "
                "dimensions ", v.dims);

            ++i;
            while (i < ids.size() && i0 - ids.safe[ids.size()-1-i] <= 1) {
                i0 = ids.safe[ids.size()-1-i];

                phypp_check(i0 != i1, "remove indices contain duplicates");
                phypp_check(i0 < osize, "trying to erase value ", i0, " in vector of "
                    "dimensions ", v.dims);

                ++i;
            }

            v.data.erase(v.data.begin()+i0*pitch, v.data.begin()+(i1+1)*pitch);
        }

        v.dims[0] -= ids.size();
    }

    template<std::size_t Dim, typename Type,
        typename enable = typename std::enable_if<!std::is_pointer<Type>::value>::type>
    vec<Dim,Type> remove(vec<Dim,Type> v, const vec1u& ids) {
        inplace_remove(v, ids);
        return v;
    }

    template<std::size_t N, std::size_t Dim, typename Type1, typename Type2 = Type1,
        typename enable = typename std::enable_if<(N < Dim) && !std::is_pointer<Type1>::value>::type>
    void append(vec<Dim,Type1>& t1, const vec<Dim,Type2>& t2) {
        if (t1.empty()) {
            t1 = t2;
            return;
        }

        if (t2.empty()) return;

        std::size_t n1 = t1.dims[N], n2 = t2.dims[N];
        phypp_check(t1.size()/n1 == t2.size()/n2, "cannot append dimension ", N, " in (", t1.dims,
            ") and (", t2.dims, ")");

        // TODO: (optimization) no need for this copy?
        auto tmp = t1;
        t1.dims[N] += n2;
        t1.resize();

        t1(repeat<N>(_), uindgen(n1), repeat<Dim-N-1>(_)) = tmp;
        t1(repeat<N>(_), n1+uindgen(n2), repeat<Dim-N-1>(_)) = t2;
    }

    template<std::size_t N, std::size_t Dim, typename Type1, typename Type2 = Type1,
        typename enable = typename std::enable_if<(N < Dim) && !std::is_pointer<Type1>::value>::type>
    void prepend(vec<Dim,Type1>& t1, const vec<Dim,Type2>& t2) {
        if (t1.empty()) {
            t1 = t2;
            return;
        }

        if (t2.empty()) return;

        std::size_t n1 = t1.dims[N], n2 = t2.dims[N];
        phypp_check(t1.size()/n1 == t2.size()/n2, "cannot prepend dimension ", N, " in (", t1.dims,
            ") and (", t2.dims, ")");

        // TODO: (optimization) no need for this copy?
        auto tmp = t1;
        t1.dims[N] += n2;
        t1.resize();

        t1(repeat<N>(_), uindgen(n2), repeat<Dim-N-1>(_)) = t2;
        t1(repeat<N>(_), n2+uindgen(n1), repeat<Dim-N-1>(_)) = tmp;
    }

    template<typename Type1, typename Type2 = Type1,
        typename enable = typename std::enable_if<!std::is_pointer<Type1>::value>::type>
    void append(vec<1,Type1>& t1, const vec<1,Type2>& t2) {
        t1.data.insert(t1.data.end(), t2.begin(), t2.end());
        t1.dims[0] += t2.dims[0];
    }

    template<typename Type1, typename Type2 = Type1,
        typename enable = typename std::enable_if<!std::is_pointer<Type1>::value>::type>
    void prepend(vec<1,Type1>& t1, const vec<1,Type2>& t2) {
        t1.data.insert(t1.data.begin(), t2.begin(), t2.end());
        t1.dims[0] += t2.dims[0];
    }
}
