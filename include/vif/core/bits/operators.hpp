#ifndef VIF_INCLUDING_CORE_VEC_BITS
#error this file is not meant to be included separately, include "vif/core/vec.hpp" instead
#endif

namespace vif {
    ////////////////////////////////////////////
    //               Operators                //
    ////////////////////////////////////////////

    // Mathematical operators
    namespace impl {
        struct op_mul_t {};
        struct op_div_t {};
        struct op_mod_t {};
        struct op_add_t {};
        struct op_sub_t {};

        struct op_node_t {
            op_mul_t operator * (op_node_t);
            op_div_t operator / (op_node_t);
            op_mod_t operator % (op_node_t);
            op_add_t operator + (op_node_t);
            op_sub_t operator - (op_node_t);
        };

        #define OP_TYPE(op) decltype(impl::op_node_t{} op impl::op_node_t{})

        template<typename T>
        using math_bake_type = typename std::decay<typename std::remove_pointer<
            meta::data_type_t<typename std::remove_pointer<T>::type>>::type>::type;

        template<typename OP, typename T, typename U>
        struct op_res_t;

        template<typename T, typename U>
        struct op_res_t<op_mul_t, T, U> {
            using type = typename std::decay<decltype(std::declval<math_bake_type<T>>() *
                std::declval<math_bake_type<U>>())>::type;
        };

        template<typename T, typename U>
        struct op_res_t<op_div_t, T, U> {
            using type = typename std::decay<decltype(std::declval<math_bake_type<T>>() /
                std::declval<math_bake_type<U>>())>::type;
        };

        template<typename T, typename U>
        struct op_res_t<op_mod_t, T, U> {
            using type = typename std::decay<decltype(std::declval<math_bake_type<T>>() %
                std::declval<math_bake_type<U>>())>::type;
        };

        template<typename T, typename U>
        struct op_res_t<op_add_t, T, U> {
            using type = typename std::decay<decltype(std::declval<math_bake_type<T>>() +
                std::declval<math_bake_type<U>>())>::type;
        };

        template<typename T, typename U>
        struct op_res_t<op_sub_t, T, U> {
            using type = typename std::decay<decltype(std::declval<math_bake_type<T>>() -
                std::declval<math_bake_type<U>>())>::type;
        };

        template<typename T>
        const T& get_element_(const T& t, uint_t i) {
            return t;
        }

        template<std::size_t Dim, typename T>
        auto get_element_(const vec<Dim,T>& t, uint_t i) -> decltype(t.safe[i]) {
            return t.safe[i];
        }
    }

    #define VECTORIZE(op, sop) \
        template<std::size_t Dim, typename T, typename U> \
        vec<Dim,typename impl::op_res_t<OP_TYPE(op),T,U>::type> operator op (const vec<Dim,T>& v, const vec<Dim,U>& u) { \
            vif_check(v.dims == u.dims, "incompatible dimensions in operator '" #op \
                "' (", v.dims, " vs ", u.dims, ")"); \
            vec<Dim,typename impl::op_res_t<OP_TYPE(op),T,U>::type> tv; tv.dims = v.dims; tv.reserve(v.size()); \
            for (uint_t i : range(v)) { \
                tv.data.push_back(impl::get_element_(v, i) op impl::get_element_(u, i)); \
            } \
            return tv; \
        } \
        template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if< \
            meta::is_scalar<U>::value>::type> \
        vec<Dim,typename impl::op_res_t<OP_TYPE(op),T,U>::type> operator op (const vec<Dim,T>& v, const U& u) { \
            vec<Dim,typename impl::op_res_t<OP_TYPE(op),T,U>::type> tv; tv.dims = v.dims; tv.reserve(v.size()); \
            for (uint_t i : range(v)) { \
                tv.data.push_back(impl::get_element_(v, i) op u); \
            } \
            return tv; \
        } \
        template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if< \
            std::is_same<typename impl::op_res_t<OP_TYPE(op),T,U>::type, T>::value>::type> \
        vec<Dim,T> operator op (vec<Dim,T>&& v, const vec<Dim,U>& u) { \
            vif_check(v.dims == u.dims, "incompatible dimensions in operator '" #op \
                "' (", v.dims, " vs ", u.dims, ")"); \
            for (uint_t i : range(v)) { \
                v.data[i] sop impl::get_element_(u, i); \
            } \
            return std::move(v); \
        } \
        template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if< \
            std::is_same<typename impl::op_res_t<OP_TYPE(op),T,U>::type, T>::value && \
            meta::is_scalar<U>::value>::type> \
        vec<Dim,T> operator op (vec<Dim,T>&& v, const U& u) { \
            for (auto& t : v) { \
                t sop u; \
            } \
            return std::move(v); \
        } \
        template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if< \
            meta::is_scalar<U>::value>::type> \
        vec<Dim,typename impl::op_res_t<OP_TYPE(op),U,T>::type> operator op (const U& u, const vec<Dim,T>& v) { \
            vec<Dim,typename impl::op_res_t<OP_TYPE(op),T,U>::type> tv; tv.dims = v.dims; tv.reserve(v.size()); \
            for (uint_t i : range(v)) { \
                tv.data.push_back(u op impl::get_element_(v, i)); \
            } \
            return tv; \
        } \
        template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if< \
            std::is_same<typename impl::op_res_t<OP_TYPE(op),U,T>::type, T>::value>::type> \
        vec<Dim,T> operator op (const vec<Dim,U>& u, vec<Dim,T>&& v) { \
            vif_check(v.dims == u.dims, "incompatible dimensions in operator '" #op \
                "' (", v.dims, " vs ", u.dims, ")"); \
            for (uint_t i : range(v)) { \
                v.data[i] = impl::get_element_(u, i) op v.data[i]; \
            } \
            return std::move(v); \
        } \
        template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if< \
            std::is_same<typename impl::op_res_t<OP_TYPE(op),U,T>::type, T>::value && \
            meta::is_scalar<U>::value>::type> \
        vec<Dim,T> operator op (const U& u, vec<Dim,T>&& v) { \
            for (auto& t : v) { \
                t = u op t; \
            } \
            return std::move(v); \
        } \
        template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if< \
            std::is_same<typename impl::op_res_t<OP_TYPE(op),T,U>::type, T>::value && \
            !std::is_pointer<U>::value>::type> \
        vec<Dim,T> operator op (vec<Dim,T>&& v, vec<Dim,U>&& u) { \
            vif_check(v.dims == u.dims, "incompatible dimensions in operator '" #op \
                "' (", v.dims, " vs ", u.dims, ")"); \
            for (uint_t i : range(v)) { \
                v.data[i] sop u.data[i]; \
            } \
            return std::move(v); \
        }

    VECTORIZE(*, *=)
    VECTORIZE(+, +=)
    VECTORIZE(/, /=)
    VECTORIZE(%, %=)
    VECTORIZE(-, -=)

    #undef VECTORIZE

    // Logical operators
    #define VECTORIZE(op) \
        template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if<!meta::is_vec<U>::value>::type> \
        vec<Dim,bool> operator op (const vec<Dim,T>& v, const U& u) { \
            vec<Dim,bool> tv(v.dims); \
            for (uint_t i : range(v)) { \
                tv.safe[i] = (v.safe[i] op u); \
            } \
            return tv; \
        } \
        \
        template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if<!meta::is_vec<U>::value>::type> \
        vec<Dim,bool> operator op (const U& u, const vec<Dim,T>& v) { \
            vec<Dim,bool> tv(v.dims); \
            for (uint_t i : range(v)) { \
                tv.safe[i] = (u op v.safe[i]); \
            } \
            return tv; \
        } \
        \
        template<std::size_t Dim, typename T, typename U> \
        vec<Dim,bool> operator op (const vec<Dim,T>& v, const vec<Dim,U>& u) { \
            vif_check(v.dims == u.dims, "incompatible dimensions in operator '" #op \
                "' (", v.dims, " vs ", u.dims, ")"); \
            vec<Dim,bool> tv(v.dims); \
            for (uint_t i : range(v)) { \
                tv.safe[i] = (v.safe[i] op u.safe[i]); \
            } \
            return tv; \
        }

    VECTORIZE(==)
    VECTORIZE(!=)
    VECTORIZE(<)
    VECTORIZE(<=)
    VECTORIZE(>)
    VECTORIZE(>=)

    #undef VECTORIZE

    namespace meta {
        template<typename T>
        using is_bool = std::is_same<typename std::decay<typename std::remove_pointer<T>::type>::type, bool>;
    }

    #define VECTORIZE(op) \
        template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if< \
            meta::is_bool<T>::value && meta::is_bool<U>::value>::type> \
        vec<Dim,bool> operator op (const vec<Dim,T>& v1, const vec<Dim,U>& v2) { \
            vif_check(v1.dims == v2.dims, "incompatible dimensions in operator '" #op \
                "' (", v1.dims, " vs ", v2.dims, ")"); \
            vec<Dim,bool> tv = v1; \
            for (uint_t i : range(v1)) { \
                tv.safe[i] = tv.safe[i] op v2.safe[i]; \
            } \
            return tv; \
        } \
        template<std::size_t Dim, typename U, typename enable = typename std::enable_if< \
            meta::is_bool<U>::value>::type> \
        vec<Dim,bool> operator op (vec<Dim,bool>&& v1, const vec<Dim,U>& v2) { \
            vif_check(v1.dims == v2.dims, "incompatible dimensions in operator '" #op \
                "' (", v1.dims, " vs ", v2.dims, ")"); \
            for (uint_t i : range(v1)) { \
                v1.safe[i] = v1.safe[i] op v2.safe[i]; \
            } \
            return std::move(v1); \
        } \
        template<std::size_t Dim, typename T, typename enable = typename std::enable_if< \
            meta::is_bool<T>::value>::type> \
        vec<Dim,bool> operator op (const vec<Dim,T>& v1, vec<Dim,bool>&& v2) { \
            vif_check(v1.dims == v2.dims, "incompatible dimensions in operator '" #op \
                "' (", v1.dims, " vs ", v2.dims, ")"); \
            for (uint_t i : range(v2)) { \
                v2.safe[i] = v1.safe[i] op v2.safe[i]; \
            } \
            return std::move(v2); \
        } \
        template<std::size_t Dim> \
        vec<Dim,bool> operator op (vec<Dim,bool>&& v1, vec<Dim,bool>&& v2) { \
            vif_check(v1.dims == v2.dims, "incompatible dimensions in operator '" #op \
                "' (", v1.dims, " vs ", v2.dims, ")"); \
            for (uint_t i : range(v1)) { \
                v1.safe[i] = v1.safe[i] op v2.safe[i]; \
            } \
            return std::move(v1); \
        } \
        template<std::size_t Dim, typename T, typename enable = typename std::enable_if< \
            meta::is_bool<T>::value>::type> \
        vec<Dim,bool> operator op (const vec<Dim,T>& v1, bool b) { \
            vec<Dim,bool> tv = v1; \
            for (uint_t i : range(v1)) { \
                tv.safe[i] = tv.safe[i] op b; \
            } \
            return tv; \
        } \
        template<std::size_t Dim, typename T, typename enable = typename std::enable_if< \
            meta::is_bool<T>::value>::type> \
        vec<Dim,bool> operator op (bool b, const vec<Dim,T>& v2) { \
            vec<Dim,bool> tv = v2; \
            for (uint_t i : range(v2)) { \
                tv.safe[i] = b op tv.safe[i]; \
            } \
            return tv; \
        } \
        template<std::size_t Dim> \
        vec<Dim,bool> operator op (vec<Dim,bool>&& v1, bool b) { \
            for (uint_t i : range(v1)) { \
                v1.safe[i] = v1.safe[i] op b; \
            } \
            return std::move(v1); \
        } \
        template<std::size_t Dim> \
        vec<Dim,bool> operator op (bool b, vec<Dim,bool>&& v2) { \
            for (uint_t i : range(v2)) { \
                v2.safe[i] = b op v2.safe[i]; \
            } \
            return std::move(v2); \
        }

    VECTORIZE(&&)
    VECTORIZE(||)

    #undef VECTORIZE

    template<std::size_t Dim, typename T, typename enable = typename std::enable_if<
        meta::is_bool<T>::value>::type>
    vec<Dim,bool> operator ! (const vec<Dim,T>& v) {
        vec<Dim,bool> tv = v;
        for (uint_t i : range(v)) {
            tv.safe[i] = !tv.safe[i];
        }

        return tv;
    }

    template<std::size_t Dim>
    vec<Dim,bool> operator ! (vec<Dim,bool>&& v) {
        for (uint_t i : range(v)) {
            v.safe[i] = !v.safe[i];
        }

        return std::move(v);
    }

    // Print a vector into a stream.
    template<typename O, std::size_t Dim, typename Type, typename enable = typename std::enable_if<!std::is_same<Type, bool>::value>::type>
    O& operator << (O& o, const vec<Dim,Type>& v) {
        o << '{';
        for (uint_t i : range(v)) {
            if (i != 0) o << ", ";
            o << v.safe[i];
        }
        o << '}';

        return o;
    }

    template<typename O, std::size_t Dim>
    O& operator << (O& o, const vec<Dim,bool>& v) {
        o << '{';
        for (uint_t i : range(v)) {
            if (i != 0) o << ", ";
            o << bool(v.safe[i]);
        }
        o << '}';

        return o;
    }
}
