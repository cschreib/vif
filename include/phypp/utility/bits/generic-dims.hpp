#ifndef PHYPP_INCLUDING_GENERIC_BITS
#error this file is not meant to be included separately, include "phypp/utilty/generic.hpp" instead
#endif

namespace phypp {
    template<typename T, typename enable = typename std::enable_if<!meta::is_vec<T>::value>::type>
    T flatten(T&& t) {
        return t;
    }

    template<std::size_t Dim, typename Type>
    vec<1,Type> flatten(const vec<Dim,Type>& v) {
        vec<1,Type> r;
        r.dims[0] = v.data.size();
        r.data = v.data;
        return r;
    }

    template<std::size_t Dim, typename Type>
    vec<1,Type> flatten(vec<Dim,Type>&& v) {
        vec<1,Type> r;
        r.dims[0] = v.data.size();
        r.data = std::move(v.data);
        return r;
    }

    template<std::size_t Dim, typename Type>
    vec<1,Type*> flatten(const vec<Dim,Type*>& v) {
        vec<1,Type*> r(impl::vec_ref_tag, v.parent);
        r.dims[0] = v.data.size();
        r.data = v.data;
        return r;
    }

    template<std::size_t Dim, typename Type>
    vec<1,Type*> flatten(vec<Dim,Type*>&& v) {
        vec<1,Type*> r(impl::vec_ref_tag, v.parent);
        r.dims[0] = v.data.size();
        r.data = std::move(v.data);
        return r;
    }

    template<std::size_t Dim, typename Type, typename ... Args>
    vec<meta::dim_total<Args...>::value, Type> reform(const vec<Dim,Type>& v, Args&& ... args) {
        vec<meta::dim_total<Args...>::value, Type> r;
        impl::set_array(r.dims, std::forward<Args>(args)...);
        uint_t nsize = 1;
        for (uint_t i : range(meta::dim_total<Args...>::value)) {
            nsize *= r.dims[i];
        }

        phypp_check(v.size() == nsize,
            "incompatible dimensions (", v.dims, " vs ", r.dims, ")");

        r.data = v.data;

        return r;
    }

    template<std::size_t Dim, typename Type, typename ... Args>
    vec<meta::dim_total<Args...>::value, Type> reform(vec<Dim,Type>&& v, Args&& ... args) {
        vec<meta::dim_total<Args...>::value, Type> r;
        impl::set_array(r.dims, std::forward<Args>(args)...);
        uint_t nsize = 1;
        for (uint_t i : range(meta::dim_total<Args...>::value)) {
            nsize *= r.dims[i];
        }

        phypp_check(v.size() == nsize,
            "incompatible dimensions ("+strn(v.dims)+" vs "+strn(r.dims)+")");

        r.data = std::move(v.data);

        return r;
    }

    template<std::size_t Dim, typename Type, typename ... Args>
    vec<meta::dim_total<Args...>::value, Type*> reform(const vec<Dim,Type*>& v, Args&& ... args) {
        vec<meta::dim_total<Args...>::value, Type*> r(impl::vec_ref_tag, v.parent);
        impl::set_array(r.dims, std::forward<Args>(args)...);
        uint_t nsize = 1;
        for (uint_t i : range(meta::dim_total<Args...>::value)) {
            nsize *= r.dims[i];
        }

        phypp_check(v.size() == nsize,
            "incompatible dimensions ("+strn(v.dims)+" vs "+strn(r.dims)+")");

        r.data = v.data;

        return r;
    }

    template<std::size_t Dim, typename Type, typename ... Args>
    vec<meta::dim_total<Args...>::value, Type*> reform(vec<Dim,Type*>&& v, Args&& ... args) {
        vec<meta::dim_total<Args...>::value, Type*> r(impl::vec_ref_tag, v.parent);
        impl::set_array(r.dims, std::forward<Args>(args)...);
        uint_t nsize = 1;
        for (uint_t i : range(meta::dim_total<Args...>::value)) {
            nsize *= r.dims[i];
        }

        phypp_check(v.size() == nsize,
            "incompatible dimensions ("+strn(v.dims)+" vs "+strn(r.dims)+")");

        r.data = std::move(v.data);

        return r;
    }
}
