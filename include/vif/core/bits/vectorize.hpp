#ifndef PHYPP_INCLUDING_CORE_VEC_BITS
#error this file is not meant to be included separately, include "phypp/core/vec.hpp" instead
#endif

namespace phypp {
    ////////////////////////////////////////////
    //         Vectorization helpers          //
    ////////////////////////////////////////////

    #define PHYPP_VECTORIZE(name) \
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

    #define PHYPP_VECTORIZE2(name) \
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

    #define PHYPP_VECTORIZE_REN(name, orig) \
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

    // Create an overloaded lambda that supports both scalars and vectors for the first argument.
    namespace impl {
        template<typename L>
        struct vectorized_lambda_first_t {
            L lambda;

            vectorized_lambda_first_t(L tlam) : lambda(tlam) {}

            template<typename T, typename ... Args, typename enable =
                typename std::enable_if<!meta::is_vec<typename std::decay<T>::type>::value>::type>
            auto operator()(T&& t, Args&& ... args) ->
                decltype(lambda(std::forward<T>(t), std::forward<Args>(args)...)) {
                return lambda(std::forward<T>(t), std::forward<Args>(args)...);
            }

            template<typename T, std::size_t D, typename ... Args>
            auto operator()(const vec<D,T>& t, Args&& ... args) ->
                vec<D,typename std::decay<decltype(lambda(std::declval<const meta::rtype_t<T>&>(),
                    std::forward<Args>(args)...))>::type> {
                vec<D,typename std::decay<decltype(lambda(std::declval<const meta::rtype_t<T>&>(),
                    std::forward<Args>(args)...))>::type> ret(t.dims);
                for (uint_t i : range(t)) {
                    ret.safe[i] = lambda(t.safe[i], std::forward<Args>(args)...);
                }
                return ret;
            }

            template<typename T, std::size_t D, typename ... Args>
            auto operator()(vec<D,T>&& t, Args&& ... args) ->
                vec<D,typename std::decay<decltype(lambda(std::declval<meta::rtype_t<T>>(),
                    std::forward<Args>(args)...))>::type> {
                vec<D,typename std::decay<decltype(lambda(std::declval<meta::rtype_t<T>>(),
                    std::forward<Args>(args)...))>::type> ret(t.dims);
                for (uint_t i : range(t)) {
                    ret.safe[i] = lambda(std::move(t.safe[i]), std::forward<Args>(args)...);
                }
                return ret;
            }
        };
    }

    template<typename T>
    impl::vectorized_lambda_first_t<typename std::decay<T>::type> vectorize_lambda_first(T&& t) {
        return impl::vectorized_lambda_first_t<typename std::decay<T>::type>(std::move(t));
    }

    // Create an overloaded lambda that supports either scalars or vectors for all arguments.
    // If multiple vectors are found in the argument list, they are iterated jointly and
    // therefore must have the same dimensions.
    namespace impl {
        template<typename T>
        struct elem_dim : std::integral_constant<std::size_t, 0> {};

        template<std::size_t D, typename T>
        struct elem_dim<vec<D,T>> : std::integral_constant<std::size_t, D> {};

        template<typename ... Args>
        using common_dim = meta::max<std::size_t, elem_dim<typename std::decay<Args>::type>::value...>;

        template <std::size_t D>
        struct has_right_dim {
            template <typename T>
            struct type : meta::bool_constant<elem_dim<T>::value == D || elem_dim<T>::value == 0> {};
        };

        template <typename T, typename U>
        void lambda_check_dims(T& dims, bool& set, const U& u) {}

        template <typename T, std::size_t D, typename U>
        void lambda_check_dims(T& dims, bool& set, const vec<D,U>& u) {
            if (!set) {
                dims = u.dims;
                set = true;
            } else {
                phypp_check(dims == u.dims, "incompatible dimensions in lambda call (",
                    dims, " vs ", u.dims, ")");
            }
        }

        template<typename L>
        struct vectorize_lambda_t {
            L lambda;

            vectorize_lambda_t(L tlam) : lambda(tlam) {}

            // Get 'i'th element from argument
            template <typename T, typename ... Args>
            static T& get(uint_t i, T& t) {
                return t;
            }

            template <typename T, typename ... Args>
            static const T& get(uint_t i, const T& t) {
                return t;
            }

            template <std::size_t D, typename T, typename ... Args>
            static auto get(uint_t i, vec<D,T>& t) -> decltype(t.safe[i]) {
                return t.safe[i];
            }

            template <std::size_t D, typename T, typename ... Args>
            static auto get(uint_t i, const vec<D,T>& t) -> decltype(t.safe[i]) {
                return t.safe[i];
            }

            template <typename ... Args>
            static void swallow(Args&&...) {}

            // Bake return type
            template <typename ... Args>
            using scalar_return_type = decltype(std::declval<L>()(get(0, std::declval<Args>())...));

            template <typename ... Args>
            using vector_return_type = vec<common_dim<Args...>::value, scalar_return_type<Args...>>;

            // Full scalar call
            // ----------------
            template <typename ... Args>
            void run(std::true_type, std::true_type, Args&& ... args) {
                lambda(std::forward<Args>(args)...);
            }

            template <typename ... Args>
            scalar_return_type<Args...> run(std::true_type, std::false_type, Args&& ... args) {
                return lambda(std::forward<Args>(args)...);
            }

            // Vectorized call
            // ----------------

            template <typename ... Args>
            void run(std::false_type, std::true_type, Args&& ... args) {
                constexpr const std::size_t D = common_dim<Args...>::value;
                static_assert(meta::are_all_true<meta::bool_list<has_right_dim<D>::template type<Args>::value...>>::value,
                    "incompatible number of dimensions in lambda call");

                bool set = false;
                std::array<std::size_t,D> dims; // only used to get the dimensions
                swallow((lambda_check_dims(dims, set, args), 0)...);

                std::size_t size = 1;
                for (uint_t i : range(D)) {
                    size *= dims[i];
                }

                for (uint_t i : range(size)) {
                    lambda(get(i, args)...);
                }
            }

            template <typename ... Args>
            vector_return_type<Args...> run(std::false_type, std::false_type, Args&& ... args) {
                vector_return_type<Args...> ret;
                constexpr const std::size_t D = common_dim<Args...>::value;
                static_assert(meta::are_all_true<meta::bool_list<has_right_dim<D>::template type<Args>::value...>>::value,
                    "incompatible number of dimensions in lambda call");

                bool set = false;
                swallow((lambda_check_dims(ret.dims, set, args), 0)...);

                ret.resize();

                for (uint_t i : range(ret)) {
                    ret.safe[i] = lambda(get(i, args)...);
                }

                return ret;
            }

            // Generic call dispatcher
            // -----------------------

            // Needs to check if:
            // 1) all the parameters are scalars, to return a scalar or a vector
            // 2) the return value is void, to avoid creating a return value
            template<typename ... Args>
            auto operator()(Args&& ... args) ->
                decltype(this->run(meta::bool_constant<common_dim<Args...>::value == 0>{},
                             std::is_same<scalar_return_type<Args...>, void>{},
                             std::forward<Args>(args)...)) {
                return this->run(meta::bool_constant<common_dim<Args...>::value == 0>{},
                           std::is_same<scalar_return_type<Args...>, void>{},
                           std::forward<Args>(args)...);
            }
        };
    }

    template<typename T>
    impl::vectorize_lambda_t<typename std::decay<T>::type> vectorize_lambda(T&& t) {
        return impl::vectorize_lambda_t<typename std::decay<T>::type>(std::move(t));
    }
}
