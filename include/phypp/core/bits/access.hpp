#ifndef PHYPP_INCLUDING_CORE_VEC_BITS
#error this file is not meant to be included separately, include "phypp/core/vec.hpp" instead
#endif

namespace phypp {
    namespace impl {
       // Tag a value to indicate that is has to repeated N times
        template<std::size_t N, typename T>
        struct repeated_value {
            T value;
        };
    }

    template<std::size_t N, typename T>
    impl::repeated_value<N,T> repeat(T t) {
        return impl::repeated_value<N,T>{t};
    }

namespace impl {
namespace vec_access {
    // Unwrap repeated_value arguments into individual arguments
    template<typename ... IArgs>
    struct unfolder;

    template<typename F, typename ... Args>
    struct deduce_return_type {
        using type = decltype(std::declval<typename std::decay<F>::type>()(std::declval<Args&>()...));
    };

    template<>
    struct unfolder<> {
        template<typename F, typename ... OArgs>
        using out_type = typename deduce_return_type<F,OArgs...>::type;

        template<typename ... OArgs>
        struct out {
            template<typename F>
            static out_type<F,OArgs...> unfold(F&& func, OArgs& ... oargs) {
                return func(oargs...);
            }
        };
    };

    template<typename T, typename ... IArgs>
    struct unfolder<T, IArgs...> {
        template<typename F, typename ... OArgs>
        using out_type = typename unfolder<IArgs...>::template out_type<F, OArgs..., T>;

        template<typename ... OArgs>
        struct out {
            template<typename F>
            static out_type<F,OArgs...> unfold(F&& func, T& t, IArgs& ... iargs, OArgs& ... oargs) {
                using actor = typename unfolder<IArgs...>::template out<OArgs..., T>;
                return actor::unfold(std::forward<F>(func), iargs..., oargs..., t);
            }
        };
    };

    template<std::size_t N, typename T, typename ... IArgs>
    struct unfolder<repeated_value<N,T>, IArgs...> {
        template<typename F, std::size_t I, typename ... OArgs>
        struct out_type_impl;

        template<typename F, typename ... OArgs>
        struct out_type_impl<F, 0, OArgs...> {
            using type = typename unfolder<IArgs...>::template out_type<F, OArgs...>;
        };

        template<typename F, std::size_t I, typename ... OArgs>
        struct out_type_impl {
            using type = typename out_type_impl<F, I-1, OArgs..., T>::type;
        };

        template<typename F, typename ... OArgs>
        using out_type = typename out_type_impl<F, N, OArgs...>::type;

        template<typename ... OArgs>
        struct out {
            template<typename F, typename U>
            static out_type<F,OArgs...> unfold_repeat(F&& func, U& t, meta::cte_t<0>, IArgs& ... iargs,
                OArgs& ... oargs) {

                using actor = typename unfolder<IArgs...>::template out<OArgs...>;
                return actor::unfold(std::forward<F>(func), iargs..., oargs...);
            }

            template<typename F, std::size_t I, typename U>
            static out_type<F,OArgs...> unfold_repeat(F&& func, U& t, meta::cte_t<I>, IArgs& ... iargs,
                OArgs& ... oargs) {

                using actor = typename unfolder::template out<OArgs..., T>;
                return actor::unfold_repeat(std::forward<F>(func), t, meta::cte_t<I-1>(),
                    iargs..., oargs..., t.value);
            }

            template<typename F, typename U>
            static out_type<F,OArgs...> unfold(F&& func, U& t, IArgs& ... iargs, OArgs& ... oargs) {
                return unfold_repeat(std::forward<F>(func), t, meta::cte_t<N>(), iargs..., oargs...);
            }
        };
    };

    template<std::size_t N, typename T, typename ... IArgs>
    struct unfolder<const repeated_value<N,T>, IArgs...> {
        template<typename F, std::size_t I, typename ... OArgs>
        struct out_type_impl;

        template<typename F, typename ... OArgs>
        struct out_type_impl<F, 0, OArgs...> {
            using type = typename unfolder<IArgs...>::template out_type<F, OArgs...>;
        };

        template<typename F, std::size_t I, typename ... OArgs>
        struct out_type_impl {
            using type = typename out_type_impl<F, I-1, OArgs..., const T>::type;
        };

        template<typename F, typename ... OArgs>
        using out_type = typename out_type_impl<F, N, OArgs...>::type;

        template<typename ... OArgs>
        struct out {
            template<typename F, typename U>
            static out_type<F,OArgs...> unfold_repeat(F&& func, U& t, meta::cte_t<0>, IArgs& ... iargs,
                OArgs& ... oargs) {

                using actor = typename unfolder<IArgs...>::template out<OArgs...>;
                return actor::unfold(std::forward<F>(func), iargs..., oargs...);
            }

            template<typename F, std::size_t I, typename U>
            static out_type<F,OArgs...> unfold_repeat(F&& func, U& t, meta::cte_t<I>, IArgs& ... iargs,
                OArgs& ... oargs) {

                using actor = typename unfolder::template out<OArgs..., const T>;
                return actor::unfold_repeat(std::forward<F>(func), t, meta::cte_t<I-1>(),
                    iargs..., oargs..., t.value);
            }

            template<typename F, typename U>
            static out_type<F,OArgs...> unfold(F&& func, U& t, IArgs& ... iargs, OArgs& ... oargs) {
                return unfold_repeat(std::forward<F>(func), t, meta::cte_t<N>(), iargs..., oargs...);
            }
        };
    };

    template<typename F, typename ... Args>
    typename unfolder<Args...>::template out_type<F> unfold(F&& func, Args& ... args) {
        using actor = typename unfolder<Args...>::template out<>;
        return actor::unfold(std::forward<F>(func), args...);
    }

    template<std::size_t Dim, typename Type>
    const vec<Dim,Type>& get_parent(const vec<Dim,Type>& v) {
        return v;
    }

    template<std::size_t Dim, typename Type>
    void* get_parent(const vec<Dim,Type*>& v) {
        return v.parent;
    }

    // Generated vector dimension given index type
    template<typename T>
    struct output_dim_ :
        std::integral_constant<std::size_t, 0> {};

    template<std::size_t N, typename T>
    struct output_dim_<vec<N,T>> :
        std::integral_constant<std::size_t, N> {};

    template<>
    struct output_dim_<impl::range_impl::full_range_t> :
        std::integral_constant<std::size_t, 1> {};
    template<>
    struct output_dim_<impl::range_impl::left_range_t> :
        std::integral_constant<std::size_t, 1> {};
    template<>
    struct output_dim_<impl::range_impl::right_range_t> :
        std::integral_constant<std::size_t, 1> {};
    template<>
    struct output_dim_<impl::range_impl::left_right_range_t> :
        std::integral_constant<std::size_t, 1> {};

    template<std::size_t N, typename T>
    struct output_dim_<impl::repeated_value<N,T>> :
        std::integral_constant<std::size_t, N*output_dim_<T>::value> {};

    template<typename T>
    using output_dim = output_dim_<typename std::decay<T>::type>;

    template<typename T, typename ... Args>
    struct result_dim : std::integral_constant<std::size_t,
        output_dim<T>::value + result_dim<Args...>::value> {};

    template<typename T>
    struct result_dim<T> : std::integral_constant<std::size_t,
        output_dim<T>::value> {};

    // Needed vector dimension given index type
    template<typename T>
    struct input_dim_ :
        std::integral_constant<std::size_t, 1> {};

    template<std::size_t N, typename T>
    struct input_dim_<vec<N,T>> :
        std::integral_constant<std::size_t, N> {};

    template<>
    struct input_dim_<impl::range_impl::full_range_t> :
        std::integral_constant<std::size_t, 1> {};
    template<>
    struct input_dim_<impl::range_impl::left_range_t> :
        std::integral_constant<std::size_t, 1> {};
    template<>
    struct input_dim_<impl::range_impl::right_range_t> :
        std::integral_constant<std::size_t, 1> {};
    template<>
    struct input_dim_<impl::range_impl::left_right_range_t> :
        std::integral_constant<std::size_t, 1> {};

    template<std::size_t N, typename T>
    struct input_dim_<impl::repeated_value<N,T>> :
        std::integral_constant<std::size_t, N*input_dim_<T>::value> {};

    template<typename T>
    using input_dim = input_dim_<typename std::decay<T>::type>;

    template<typename T, typename ... Args>
    struct accessed_dim : std::integral_constant<std::size_t,
        input_dim<T>::value + accessed_dim<Args...>::value> {};

    template<typename T>
    struct accessed_dim<T> : std::integral_constant<std::size_t,
        input_dim<T>::value> {};

    template<typename T>
    struct is_index_vector : std::false_type {};

    template<std::size_t Dim, typename T>
    struct is_index_vector<vec<Dim,T>> : std::integral_constant<bool,
        std::is_integral<typename std::decay<typename std::remove_pointer<T>::type>::type>::value> {};

    template<typename T>
    struct is_index_base : std::integral_constant<bool,
        std::is_integral<T>::value || meta::is_range<T>::value || is_index_vector<T>::value> {};

    template<typename T>
    struct is_repeated_index : std::false_type {};

    template<std::size_t N, typename T>
    struct is_repeated_index<impl::repeated_value<N,T>> : is_index_base<T> {};

    template<typename T>
    struct is_index : std::integral_constant<bool,
        is_index_base<T>::value || is_repeated_index<T>::value> {};

    template<typename ... Args>
    struct are_indices : std::integral_constant<bool,
        meta::are_all_true<meta::bool_list<is_index<Args>::value...>>::value> {};

    template<>
    struct are_indices<> : std::true_type {};

    inline impl::range_impl::range_t<uint_t> range(impl::range_impl::full_range_t, uint_t size) {
        return phypp::range(size);
    }

    inline impl::range_impl::range_t<uint_t> range(const impl::range_impl::left_range_t& rng, uint_t size) {
        phypp::impl::range_impl::check_bounds(rng, size);
        return phypp::range(rng.last+1);
    }

    inline impl::range_impl::range_t<uint_t> range(const impl::range_impl::right_range_t& rng, uint_t size) {
        phypp::impl::range_impl::check_bounds(rng, size);
        return phypp::range(rng.first, size);
    }

    inline impl::range_impl::range_t<uint_t> range(const impl::range_impl::left_right_range_t& rng, uint_t size) {
        phypp::impl::range_impl::check_bounds(rng, size);
        return phypp::range(rng.first, rng.last+1);
    }

    template<typename T>
    uint_t to_idx_(uint_t size, T ui, meta::cte_t<false>) {
        phypp_check(ui < size, "operator[]: index out of bounds (", ui, " vs. ", size, ")");
        return ui;
    }

    template<typename T>
    uint_t to_idx_(uint_t size, T i, meta::cte_t<true>) {
        while (i < 0) i += size;
        uint_t ui(i);
        phypp_check(ui < size, "operator[]: index out of bounds (", ui, " vs. ", size, ")");
        return ui;
    }

    template<std::size_t I, std::size_t D, typename T>
    uint_t to_idx_(std::array<uint_t,D> dims, T ui, meta::cte_t<false>) {
        phypp_check(ui < dims[I], "operator(): index out of bounds (", ui, " vs. ", dims[D], ")");
        return ui;
    }

    template<std::size_t I, std::size_t D, typename T>
    uint_t to_idx_(std::array<uint_t,D> dims, T i, meta::cte_t<true>) {
        while (i < 0) i += dims[I];
        uint_t ui(i);
        phypp_check(ui < dims[I], "operator(): index out of bounds (", ui, " vs. ", dims[I], ")");
        return ui;
    }

    template<typename T, typename enable =
        typename std::enable_if<std::is_integral<T>::value>::type>
    uint_t to_idx(uint_t size, T ix) {
        return to_idx_(size, ix, meta::cte_t<std::is_signed<T>::value>());
    }

    template<std::size_t I, std::size_t D, typename T, typename enable =
        typename std::enable_if<std::is_integral<T>::value>::type>
    uint_t to_idx(std::array<uint_t,D> dims, T ix) {
        return to_idx_<I>(dims, ix, meta::cte_t<std::is_signed<T>::value>());
    }

    template<std::size_t D>
    uint_t pitch(std::array<uint_t,D> dims, uint_t i) {
        uint_t p = 1;
        for (uint_t j = i+1; j < D; ++j) {
            p *= dims[j];
        }
        return p;
    }

    // Helper to build the result of v(_, ids, 5), i.e. when at least one index is not scalar.
    // The result is another array.
    template<bool IsSafe, bool IsConst, std::size_t Dim, std::size_t ODim, typename Type,
        typename ... Args>
    struct helper_ {
        static_assert(are_indices<Args...>::value, "vector access can only be made with "
            "integers or '_'");

        using rptype = typename std::remove_pointer<Type>::type;
        using itype = typename std::conditional<IsConst,
            const vec<Dim,Type>, vec<Dim,Type>>::type;
        using type = typename std::conditional<IsConst,
            vec<ODim, const rptype*>, vec<ODim, rptype*>>::type;

        // Functions to build the dimension of the resulting vector
        template<std::size_t IT, std::size_t IV, typename T>
        static void do_resize_impl_(type&, itype&, meta::cte_t<IT>, meta::cte_t<IV>, const T&, std::false_type) {}

        template<std::size_t IT, std::size_t IV, typename T>
        static void do_resize_impl_(type& t, itype& v, meta::cte_t<IT>, meta::cte_t<IV>, const T& rng, std::true_type) {
            t.dims[IT] = impl::range_impl::range_size(rng, v.dims[IV]);
        }

        template<std::size_t IT, std::size_t IV, typename T>
        static void do_resize_(type& t, itype& v, meta::cte_t<IT> d1, meta::cte_t<IV> d2, const T& i) {
            do_resize_impl_(t, v, d1, d2, i, meta::is_range<T>{});
        }

        template<std::size_t IT, std::size_t IV, typename T>
        static void do_resize_(type& t, itype&, meta::cte_t<IT>, meta::cte_t<IV>, const vec<1,T>& ids) {
            t.dims[IT] = ids.size();
        }

        template<std::size_t IT, std::size_t IV, typename T, typename ... Args2>
        static void resize_(type& t, itype& v, meta::cte_t<IT> it, meta::cte_t<IV> iv, const T& i,
            const Args2& ... args) {
            do_resize_(t, v, it, iv, i);
            resize_(t, v, meta::cte_t<IT+output_dim<T>::value>(), meta::cte_t<IV+input_dim<T>::value>(), args...);
        }

        static void resize_(type& t, itype& v, meta::cte_t<ODim>, meta::cte_t<Dim>) {
            t.resize();
        }

        // Adapter to switch between safe/unsafe array indexing
        template<std::size_t D, typename T>
        static uint_t to_idx_(itype& v, const T& t, std::true_type) {
            return v.template to_idx<D>(t);
        }

        template<std::size_t D>
        static uint_t to_idx_(itype&, uint_t t, std::false_type) {
            return t;
        }

        template<std::size_t D, typename T>
        static uint_t to_idx(itype& v, const T& t) {
            return to_idx_<D>(v, t, std::integral_constant<bool,IsSafe>{});
        }

        // Functions to populate the resulting vector
        template<typename T>
        static void make_indices_impl_(type& t, uint_t& itx, itype& v, uint_t ivx, meta::cte_t<Dim-1>,
            const std::array<uint_t, Dim>& pitch, std::false_type, const T& ix) {

            t.data[itx] = impl::ptr<Type>(v.data[ivx+to_idx<Dim-1>(v,ix)]);
            ++itx;
        }

        template<typename T>
        static void make_indices_impl_(type& t, uint_t& itx, itype& v, uint_t ivx, meta::cte_t<Dim-1>,
            const std::array<uint_t, Dim>& pitch, std::true_type, const T& rng) {

            for (uint_t j : range(rng, v.dims[Dim-1])) {
                t.data[itx] = impl::ptr<Type>(v.data[ivx+j]);
                ++itx;
            }
        }

        template<typename T>
        static void make_indices_impl_(type& t, uint_t& itx, itype& v, uint_t ivx, meta::cte_t<Dim-1> d,
            const std::array<uint_t, Dim>& pitch, const T& ix) {

            make_indices_(t, itx, v, ivx, d, pitch, meta::is_range<T>{}, ix);
        }

        template<typename T>
        static void make_indices_(type& t, uint_t& itx, itype& v, uint_t ivx, meta::cte_t<Dim-1>,
            const std::array<uint_t, Dim>& pitch, const vec<1,T>& ids) {

            for (uint_t j : ids) {
                t.data[itx] = impl::ptr<Type>(v.data[ivx+to_idx<Dim-1>(v,j)]);
                ++itx;
            }
        }

        template<std::size_t IV, typename T, typename ... Args2>
        static void make_indices_impl_(type& t, uint_t& itx, itype& v, uint_t ivx, meta::cte_t<IV>,
            const std::array<uint_t, Dim>& pitch, std::false_type, const T& ix, const Args2& ... i) {

            make_indices_(t, itx, v, ivx +
                to_idx<IV>(v,ix)*pitch[IV], meta::cte_t<IV+1>(), pitch, i...
            );
        }

        template<std::size_t IV, typename T, typename ... Args2>
        static void make_indices_impl_(type& t, uint_t& itx, itype& v, uint_t ivx, meta::cte_t<IV>,
            const std::array<uint_t, Dim>& pitch, std::true_type, const T& rng, const Args2& ... i) {

            for (uint_t j : range(rng, v.dims[IV])) {
                make_indices_(t, itx, v, ivx + j*pitch[IV], meta::cte_t<IV+1>(), pitch, i...);
            }
        }

        template<std::size_t IV, typename T, typename ... Args2>
        static void make_indices_(type& t, uint_t& itx, itype& v, uint_t ivx, meta::cte_t<IV> d,
            const std::array<uint_t, Dim>& pitch, const T& ix, const Args2& ... i) {

            make_indices_impl_(t, itx, v, ivx, d, pitch, meta::is_range<T>{}, ix, i...);
        }

        template<std::size_t IV, typename T, typename ... Args2>
        static void make_indices_(type& t, uint_t& itx, itype& v, uint_t ivx, meta::cte_t<IV>,
            const std::array<uint_t, Dim>& pitch, const vec<1,T>& ids, const Args2& ... i) {

            for (uint_t j : ids) {
                make_indices_(t, itx, v, ivx + to_idx<IV>(v,j)*pitch[IV], meta::cte_t<IV+1>(),
                    pitch, i...
                );
            }
        }

        template<typename ... UArgs>
        static type access_(itype& v, const UArgs& ... i) {
            type t(impl::vec_ref_tag, get_parent(v));
            resize_(t, v, meta::cte_t<0>(), meta::cte_t<0>(), i...);

            // TODO: (optimization) cache this on vector construction
            std::array<uint_t, Dim> pitch;
            for (uint_t j = 0; j < Dim; ++j) {
                pitch[j] = 1;
                for (uint_t k = j+1; k < Dim; ++k) {
                    pitch[j] *= v.dims[k];
                }
            }

            uint_t itx = 0;
            make_indices_(t, itx, v, 0, meta::cte_t<0>(), pitch, i...);
            return t;
        }

        template<typename ... UArgs>
        type operator() (itype& v, const UArgs& ... i) const {
            return access_(v, i...);
        }

        static type access(itype& v, const Args& ... i) {
            return unfold(helper_(), v, i...);
        }
    };

    template<bool IsSafe, bool IsConst, typename Type, typename Arg>
    struct helper_<IsSafe, IsConst, 1, 1, Type, Arg> {
        static_assert(is_index<Arg>::value, "vector access can only be made with "
            "integers or '_'");

        using rptype = typename std::remove_pointer<Type>::type;
        using itype = typename std::conditional<IsConst,
            const vec<1,Type>, vec<1,Type>>::type;
        using type = typename std::conditional<IsConst,
            vec<1, const rptype*>, vec<1, rptype*>>::type;

        // Adapter to switch between safe/unsafe array indexing
        template<typename T>
        static uint_t to_idx_(itype& v, const T& t, std::true_type) {
            return v.to_idx(t);
        }

        static uint_t to_idx_(itype&, uint_t t, std::false_type) {
            return t;
        }

        template<typename T>
        static uint_t to_idx(itype& v, const T& t) {
            return to_idx_(v, t, std::integral_constant<bool,IsSafe>{});
        }

        // Access functions
        template<typename T, typename enable = typename std::enable_if<meta::is_range<T>::value>::type>
        static type access(itype& v, const T& rng) {
            type t(impl::vec_ref_tag, get_parent(v));
            t.dims[0] = impl::range_impl::range_size(rng, v.dims[0]);
            t.resize();

            uint_t itx = 0;
            for (uint_t i : range(rng, v.dims[0])) {
                t.data[itx] = impl::ptr<Type>(v.data[i]);
                ++itx;
            }

            return t;
        }

        template<typename T>
        static type access(itype& v, const vec<1,T>& ids) {
            type t(impl::vec_ref_tag, get_parent(v));
            t.dims[0] = ids.size();
            t.resize();

            uint_t itx = 0;
            for (uint_t i : ids) {
                t.data[itx] = impl::ptr<Type>(v.data[to_idx(v,i)]);
                ++itx;
            }

            return t;
        }
    };

    // Helper to build the result of v(1,2,3), i.e. when all indices are scalars.
    // The result is a scalar.
    template<bool IsSafe, bool IsConst, std::size_t Dim, typename Type, typename ... Args>
    struct helper_<IsSafe, IsConst, Dim, 0, Type, Args...>  {
        static_assert(are_indices<Args...>::value, "vector access can only be made with "
            "integers or '_'");

        using rptype = typename std::remove_pointer<Type>::type;
        using itype = typename std::conditional<IsConst,
            const vec<Dim,Type>, vec<Dim,Type>>::type;
        using type = typename std::conditional<IsConst,
            const meta::dtype_t<rptype>&, meta::dtype_t<rptype>&>::type;

        // Adapter to switch between safe/unsafe array indexing
        template<std::size_t D, typename T>
        static uint_t to_idx_(itype& v, const T& t, std::true_type) {
            return v.template to_idx<D>(t);
        }

        template<std::size_t D>
        static uint_t to_idx_(itype&, uint_t t, std::false_type) {
            return t;
        }

        template<std::size_t D, typename T>
        static uint_t to_idx(itype& v, const T& t) {
            return to_idx_<D>(v, t, std::integral_constant<bool,IsSafe>{});
        }

        // Access functions
        template<typename V>
        static uint_t get_index_(V& vi, uint_t idx) {
            return idx;
        }

        template<typename V, typename T, typename ... Args2>
        static uint_t get_index_(V& v, uint_t idx, const T& ix, const Args2& ... i) {
            uint_t pitch = 1;

            // TODO: (optimization) cache this on vector construction
            for (uint_t j = Dim-sizeof...(Args2); j < Dim; ++j) {
                pitch *= v.dims[j];
            }

            idx += to_idx<Dim-1-sizeof...(Args2)>(v, ix)*pitch;

            return get_index_(v, idx, i...);
        }

        static type access(itype& v, const Args& ... i) {
            return reinterpret_cast<type>(v.safe[get_index_(v, 0, i...)]);
        }
    };

    // Helper to build the result of v(1), i.e. when all indices are scalars.
    // The result is a scalar.
    template<bool IsSafe, bool IsConst, typename Type, typename T>
    struct helper_<IsSafe, IsConst, 1, 0, Type, T> {
        static_assert(is_index<T>::value, "vector access can only be made with "
            "integers or '_'");

        using rptype = typename std::remove_pointer<Type>::type;
        using itype = typename std::conditional<IsConst,
            const vec<1,Type>, vec<1,Type>>::type;
        using type = typename std::conditional<IsConst,
            const meta::dtype_t<rptype>&, meta::dtype_t<rptype>&>::type;

        // Adapter to switch between safe/unsafe array indexing
        static uint_t to_idx_(itype& v, const T& t, std::true_type) {
            return v.to_idx(t);
        }

        static uint_t to_idx_(itype&, uint_t t, std::false_type) {
            return t;
        }

        static uint_t to_idx(itype& v, const T& t) {
            return to_idx_(v, t, std::integral_constant<bool,IsSafe>{});
        }

        // Access function
        static type access(itype& v, const T& i) {
            return v.safe[to_idx(v, i)];
        }
    };

    template<bool IsSafe, bool IsConst, std::size_t Dim, typename Type, typename ... Args>
    using helper = helper_<IsSafe, IsConst, Dim, result_dim<Args...>::value, Type, Args...>;

    template<typename V, typename T>
    auto bracket_access(V& parent, const T& rng) ->
        vec<1,meta::constify<typename V::rtype, V>*> {
        vec<1,meta::constify<typename V::rtype, V>*> v(impl::vec_ref_tag, parent);
        v.dims[0] = impl::range_impl::range_size(rng, parent.size());
        v.data.resize(v.dims[0]);

        uint_t itx = 0;
        for (uint_t i : range(rng, parent.size())) {
            v.data[itx] = impl::ptr<typename V::rtype>(parent.data[i]);
            ++itx;
        }

        return v;
    }

    template<typename ... Args>
    using count_placeholder = meta::count<meta::binary_second_apply_type_to_value_list<
        meta::type_list<typename std::decay<Args>::type...>, bool, std::is_same, placeholder_t>>;

    template<typename V, typename P, typename ... Args>
    void get_stride_offset_(V&, P&, P&, bool, uint_t) {}

    template<typename V, typename P, typename T, typename ... Args>
    void get_stride_offset_(V& parent, P& bp, P& ep, bool found, uint_t d,
        const T& id, const Args& ... i);

    template<typename V, typename P, typename ... Args>
    void get_stride_offset_(V& parent, P& bp, P& ep, bool, uint_t d,
        const placeholder_t&, const Args& ... i) {
        uint_t pitch = 1;
        for (uint_t k = d+1; k < parent.dims.size(); ++k) {
            pitch *= parent.dims[k];
        }

        bp.stride = pitch;
        ep.stride = pitch;
        ep.offset += parent.dims[d]*pitch;

        get_stride_offset_(parent, bp, ep, true, d+1, i...);
    }

    template<typename V, typename P, typename T, typename ... Args>
    void get_stride_offset_(V& parent, P& bp, P& ep, bool found, uint_t d,
        const T& id, const Args& ... i) {
        uint_t pitch = 1;
        for (uint_t k = d+1; k < parent.dims.size(); ++k) {
            pitch *= parent.dims[k];
        }

        bp.offset += id*pitch;
        ep.offset += id*pitch;
        get_stride_offset_(parent, bp, ep, found, d+1, i...);
    }

    template<typename V, typename P, typename ... Args>
    void get_stride_offset(V& parent, P& bp, P& ep, const Args& ... i) {
        get_stride_offset_(parent, bp, ep, false, 0, i...);
    }

    // Utility to allow range-based loops on strides
    template<typename T, typename P>
    struct strided_range {
        using policy = P;
        using vec = meta::decay_t<T>;
        using iterator = typename std::conditional<std::is_const<T>::value,
            typename impl::vec_iterator_type<vec,policy>::const_iterator,
            typename impl::vec_iterator_type<vec,policy>::iterator>::type;

        T& parent;
        policy begin_policy;
        policy end_policy;

        template<typename ... Args>
        strided_range(T& v, const Args& ... i) : parent(v) {
            static_assert(impl::vec_access::accessed_dim<Args...>::value == meta::vec_dim<vec>::value,
                "wrong number of indices for this vector");
            static_assert(impl::vec_access::count_placeholder<Args...>::value <= 1,
                "strides can only iterate over one dimension at a time");
            static_assert(impl::vec_access::count_placeholder<Args...>::value >= 1,
                "stride(...) must contain at least one placeholder '_'");
            impl::vec_access::get_stride_offset(parent, begin_policy, end_policy, i...);
        }

        template<typename ... Args>
        iterator begin() {
            return parent.begin(begin_policy);
        }

        template<typename ... Args>
        iterator end() {
            return parent.begin(end_policy);
        }
    };
}
}
}
