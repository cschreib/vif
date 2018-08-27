#ifndef VIF_INCLUDING_MATH_MATRIX_BITS
#error this file is not meant to be included separately, include "vif/math/matrix.hpp" instead
#endif

namespace vif {
namespace matrix {
    template<typename Type>
    struct mat;

    template<typename Type>
    struct wrapper;

    template<typename Type>
    struct const_wrapper;
}

namespace impl {
    namespace meta_impl {
        template<typename T>
        struct is_matrix : std::false_type {};

        template<typename T>
        struct is_matrix<matrix::wrapper<T>> : std::true_type {};

        template<typename T>
        struct is_matrix<matrix::const_wrapper<T>> : std::true_type {};

        template<typename T>
        struct is_matrix<matrix::mat<T>> : std::true_type {};

        template<typename T>
        struct is_matrix_vec : std::false_type {};

        template<typename T>
        struct is_matrix_vec<vec<2,T>> : std::true_type {};

        // Glue for meta::

        template<typename T>
        struct vec_dim_<matrix::wrapper<T>> : std::integral_constant<std::size_t,2> {};

        template<typename T>
        struct vec_dim_<matrix::const_wrapper<T>> : std::integral_constant<std::size_t,2> {};

        template<typename T>
        struct vec_dim_<matrix::mat<T>> : std::integral_constant<std::size_t,2> {};
    }
}

namespace meta {
    template<typename T>
    using is_matrix = impl::meta_impl::is_matrix<typename std::decay<T>::type>;

    template<typename T>
    using is_matrix_vec = impl::meta_impl::is_matrix_vec<typename std::decay<T>::type>;

    template<typename T>
    using is_matrix_like = std::integral_constant<bool, is_matrix<T>::value || is_matrix_vec<T>::value>;

    template<typename T>
    struct data_type<matrix::wrapper<T>> {
        using type = T;
    };

    template<typename T>
    struct data_type<matrix::const_wrapper<T>> {
        using type = T;
    };

    template<typename T>
    struct data_type<matrix::mat<T>> {
        using type = T;
    };

    template<typename T>
    uint_t size(const matrix::wrapper<T>& v) {
        return v.size();
    }

    template<typename T>
    uint_t size(const matrix::const_wrapper<T>& v) {
        return v.size();
    }

    template<typename T>
    uint_t size(const matrix::mat<T>& v) {
        return v.size();
    }
}

namespace matrix {
    template<typename Type>
    struct wrapper {
        using vtype = vec<2,Type>;

        vtype& base;
        std::array<uint_t,2>& dims;
        typename vtype::safe_proxy& safe;

        // Default constructor
        explicit wrapper(vtype& v) : base(v), dims(base.dims), safe(base.safe) {}

        // Copy constructor
        wrapper(const wrapper& m) = delete;

        // Move constructor
        wrapper(wrapper&& m) : base(m.base), dims(base.dims), safe(base.safe) {};

        // Import access operators
        template<typename Arg>
        auto operator [] (Arg&& arg) -> decltype(base[std::forward<Arg>(arg)]) {
            return base[std::forward<Arg>(arg)];
        }

        template<typename ... Args>
        auto operator () (Args&& ... args) -> decltype(base(std::forward<Args>(args)...)) {
            return base(std::forward<Args>(args)...);
        }

        template<typename Arg>
        auto operator [] (Arg&& arg) const -> decltype(base[std::forward<Arg>(arg)]) {
            return base[std::forward<Arg>(arg)];
        }

        template<typename ... Args>
        auto operator () (Args&& ... args) const -> decltype(base(std::forward<Args>(args)...)) {
            return base(std::forward<Args>(args)...);
        }

        // Prefix operators
        wrapper operator + () const {
            return *this;
        }

        mat<Type> operator - () const {
            mat<Type> v(base);
            for (auto& t : v) {
                t = -t;
            }

            return v;
        }

        // Forward assignment
        template<typename T>
        typename std::enable_if<meta::is_matrix<T>::value, wrapper&>::type operator = (T&& arg) {
            base = std::forward<T>(arg).base;
            return *this;
        }

        // Import usefull functions
        template<typename ... Args>
        void resize(Args&& ... args) {
            base.resize(std::forward<Args>(args)...);
        }

        uint_t size() const {
            return base.size();
        }

        bool empty() const {
            return base.empty();
        }

        typename vtype::dtype* raw_data() {
            return base.raw_data();
        }

        const typename vtype::dtype* raw_data() const {
            return base.raw_data();
        }

        typename vtype::iterator begin() {
            return base.begin();
        }

        typename vtype::const_iterator begin() const {
            return base.begin();
        }

        typename vtype::iterator end() {
            return base.end();
        }

        typename vtype::const_iterator end() const {
            return base.end();
        }

        // Implicit conversions to vector
        operator vtype& () {
            return base;
        }

        operator const vtype& () const {
            return base;
        }
    };

    template<typename Type>
    struct const_wrapper {
        using vtype = vec<2,Type>;

        const vtype& base;
        const std::array<uint_t,2>& dims;
        const typename vtype::safe_proxy& safe;

        // Default constructor
        explicit const_wrapper(const vtype& v) : base(v), dims(base.dims), safe(base.safe) {}

        // Copy constructor
        const_wrapper(const const_wrapper& m) = delete;

        // Move constructor
        const_wrapper(const_wrapper&& m) : base(m.base), dims(base.dims), safe(base.safe) {};

        // Import access operators
        template<typename Arg>
        auto operator [] (Arg&& arg) const -> decltype(base[std::forward<Arg>(arg)]) {
            return base[std::forward<Arg>(arg)];
        }

        template<typename ... Args>
        auto operator () (Args&& ... args) const -> decltype(base(std::forward<Args>(args)...)) {
            return base(std::forward<Args>(args)...);
        }

        // Prefix operators
        const_wrapper operator + () const {
            return *this;
        }

        mat<Type> operator - () const {
            mat<Type> v(base);
            for (auto& t : v) {
                t = -t;
            }

            return v;
        }

        // Disable assignment
        template<typename Arg>
        const_wrapper& operator = (Arg&& arg) = delete;

        // Import usefull functions
        uint_t size() const {
            return base.size();
        }

        bool empty() const {
            return base.empty();
        }

        const typename vtype::dtype* raw_data() const {
            return base.raw_data();
        }

        typename vtype::const_iterator begin() const {
            return base.begin();
        }

        typename vtype::const_iterator end() const {
            return base.end();
        }

        // Implicit conversion to vector
        operator const vtype& () const {
            return base;
        }
    };

    template<typename Type>
    struct mat {
        using vtype = vec<2,Type>;

        vtype base;
        std::array<uint_t,2>& dims;
        typename vtype::safe_proxy& safe;

        // Default constructor
        mat() : dims(base.dims), safe(base.safe) {}

        // Copy constructor
        mat(const mat& m) : base(m.base), dims(base.dims), safe(base.safe) {}

        // Move constructor
        mat(mat&& m) : base(std::move(m.base)), dims(base.dims), safe(base.safe) {}

        // Move in an existing vector
        explicit mat(vtype&& b) : base(std::move(b)), dims(base.dims), safe(base.safe) {}

        // Copy an existing vector
        template<typename T>
        explicit mat(const vec<2,T>& b) : base(b), dims(base.dims), safe(base.safe) {}

        // Copy an existing vector (already tagged as matrix)
        template<typename T>
        mat(const wrapper<T>& b) : base(b.base), dims(base.dims), safe(base.safe) {}

        // Copy an existing vector (already tagged as matrix)
        template<typename T>
        mat(const const_wrapper<T>& b) : base(b.base), dims(base.dims), safe(base.safe) {}

        // Dimension constructor
        template<typename ... Args, typename enable =
            typename std::enable_if<meta::is_dim_list<Args...>::value>::type>
        explicit mat(Args&& ... d) : base(std::forward<Args>(d)...), dims(base.dims), safe(base.safe) {}

        // Import access operators
        template<typename Arg>
        auto operator [] (Arg&& arg) -> decltype(base[std::forward<Arg>(arg)]) {
            return base[std::forward<Arg>(arg)];
        }

        template<typename ... Args>
        auto operator () (Args&& ... args) -> decltype(base(std::forward<Args>(args)...)) {
            return base(std::forward<Args>(args)...);
        }

        template<typename Arg>
        auto operator [] (Arg&& arg) const -> decltype(base[std::forward<Arg>(arg)]) {
            return base[std::forward<Arg>(arg)];
        }

        template<typename ... Args>
        auto operator () (Args&& ... args) const -> decltype(base(std::forward<Args>(args)...)) {
            return base(std::forward<Args>(args)...);
        }

        // Prefix operators
        mat operator + () const {
            return *this;
        }

        mat operator - () const {
            mat v(*this);
            for (auto& t : v) {
                t = -t;
            }

            return v;
        }

        // Self modifying operators
        template<typename T>
        typename std::enable_if<meta::is_matrix<T>::value, mat&>::type operator *= (T&& t) {
            *this = operator*(*this, std::forward<T>(t));
            return *this;
        }

        template<typename T>
        typename std::enable_if<meta::is_scalar<T>::value, mat&>::type operator *= (const T& t) {
            base *= t;
            return *this;
        }

        template<typename T>
        typename std::enable_if<meta::is_scalar<T>::value, mat&>::type operator /= (const T& t) {
            base /= t;
            return *this;
        }

        template<typename T>
        typename std::enable_if<meta::is_matrix<T>::value, mat&>::type operator += (const T& t) {
            base += t.base;
            return *this;
        }

        template<typename T>
        typename std::enable_if<meta::is_scalar<T>::value, mat&>::type operator += (const T& t) {
            base += t;
            return *this;
        }

        template<typename T>
        typename std::enable_if<meta::is_matrix<T>::value, mat&>::type operator -= (const T& t) {
            base -= t.base;
            return *this;
        }

        template<typename T>
        typename std::enable_if<meta::is_scalar<T>::value, mat&>::type operator -= (const T& t) {
            base -= t;
            return *this;
        }

        // Assignment
        template<typename T>
        typename std::enable_if<meta::is_matrix<T>::value, mat&>::type operator = (T&& arg) {
            base = std::forward<T>(arg).base;
            return *this;
        }

        mat& operator = (const mat& arg) {
            base = arg.base;
            return *this;
        }

        template<typename T>
        mat& operator = (mat&& arg) {
            base = std::move(arg).base;
            return *this;
        }

        // Import usefull functions
        template<typename ... Args>
        void resize(Args&& ... args) {
            base.resize(std::forward<Args>(args)...);
        }

        uint_t size() const {
            return base.size();
        }

        bool empty() const {
            return base.empty();
        }

        typename vtype::dtype* raw_data() {
            return base.raw_data();
        }

        const typename vtype::dtype* raw_data() const {
            return base.raw_data();
        }

        typename vtype::iterator begin() {
            return base.begin();
        }

        typename vtype::const_iterator begin() const {
            return base.begin();
        }

        typename vtype::iterator end() {
            return base.end();
        }

        typename vtype::const_iterator end() const {
            return base.end();
        }

        // Implicit conversions to vector
        operator vtype& () {
            return base;
        }

        operator const vtype& () const {
            return base;
        }
    };

    using mat2d = mat<double>;
    using mat2f = mat<float>;

    template<typename Type>
    mat<Type> as_matrix(vec<2,Type>&& v) {
        return mat<Type>(std::move(v));
    }

    template<typename Type>
    wrapper<Type> as_matrix(vec<2,Type>& v) {
        return wrapper<Type>(v);
    }

    template<typename Type>
    const_wrapper<Type> as_matrix(const vec<2,Type>& v) {
        return const_wrapper<Type>(v);
    }
}

namespace matrix {
    // matrix * matrix
    template<typename TypeA, typename TypeB, typename enable = typename std::enable_if<
        meta::is_matrix<TypeA>::value && meta::is_matrix<TypeB>::value
    >::type>
    auto operator * (const TypeA& a, const TypeB& b) -> mat<decltype(a(0,0)*b(0,0))> {
        vif_check(a.dims[1] == b.dims[0], "incompatible dimensions in matrix-matrix multiplication "
            "(", a.dims, " x ", b.dims, ")");

        using ntype_t = decltype(a(0,0)*b(0,0));
        mat<ntype_t> r(a.dims[0],b.dims[1]);
        for (uint_t i : range(a.dims[0]))
        for (uint_t k : range(a.dims[1]))
        for (uint_t j : range(b.dims[1])) {
            r.safe(i,j) += a.safe(i,k)*b.safe(k,j);
        }

        return r;
    }

    // matrix * 1D vector
    template<typename TypeA, typename TypeB, typename enable = typename std::enable_if<
        meta::is_matrix<TypeA>::value
    >::type>
    auto operator * (const TypeA& a, const vec<1,TypeB>& b) -> vec<1,decltype(a(0,0)*b(0,0))> {
        vif_check(a.dims[1] == b.dims[0], "incompatible dimensions in matrix-vector multiplication "
            "(", a.dims, " x ", b.dims, ")");

        using ntype_t = decltype(a(0,0)*b(0,0));
        vec<1,ntype_t> r(a.dims[0]);
        for (uint_t i : range(a.dims[0]))
        for (uint_t k : range(a.dims[1])) {
            r.safe(i) += a.safe(i,k)*b.safe(k);
        }

        return r;
    }

    // 1D vector * matrix
    template<typename TypeA, typename TypeB, typename enable = typename std::enable_if<
        meta::is_matrix<TypeA>::value
    >::type>
    auto operator * (const vec<1,TypeB>& a, const TypeA& b) -> vec<1,decltype(a(0,0)*b(0,0))> {
        vif_check(a.dims[0] == b.dims[0], "incompatible dimensions in vector-matrix multiplication "
            "(", a.dims, " x ", b.dims, ")");

        using ntype_t = decltype(a(0,0)*b(0,0));
        vec<1,ntype_t> r(b.dims[1]);
        for (uint_t k : range(b.dims[0]))
        for (uint_t i : range(b.dims[1])) {
            r.safe(i) += a.safe(k)*b.safe(k,i);
        }

        return r;
    }

    // matrix * scalar
    template<typename TypeA, typename TypeB, typename enable = typename std::enable_if<
        meta::is_matrix<TypeA>::value && std::is_arithmetic<typename std::decay<TypeB>::type>::value
    >::type>
    auto operator * (TypeA&& a, const TypeB& b) -> mat<decltype(a(0,0)*b)> {
        return mat<decltype(a(0,0)*b)>(std::forward<TypeA>(a).base*b);
    }

    // scalar * matrix
    template<typename TypeA, typename TypeB, typename enable = typename std::enable_if<
        meta::is_matrix<TypeB>::value && std::is_arithmetic<typename std::decay<TypeA>::type>::value
    >::type>
    auto operator * (const TypeB& a, TypeA&& b) -> mat<decltype(a*b(0,0))> {
        return mat<decltype(a*b(0,0))>(a*std::forward<TypeB>(b).base);
    }

    // matrix / scalar
    template<typename TypeA, typename TypeB, typename enable = typename std::enable_if<
        meta::is_matrix<TypeA>::value && std::is_arithmetic<typename std::decay<TypeB>::type>::value
    >::type>
    auto operator / (const TypeA& a, const TypeB& b) -> mat<decltype(a(0,0)/b)> {
        return mat<decltype(a(0,0)/b)>(std::forward<TypeA>(a).base/b);
    }

    // matrix + matrix
    template<typename TypeA, typename TypeB, typename enable = typename std::enable_if<
        meta::is_matrix<TypeA>::value && meta::is_matrix<TypeB>::value
    >::type>
    auto operator + (TypeA&& a, TypeB&& b) -> mat<decltype(a(0,0)+b(0,0))> {
        return mat<decltype(a(0,0)+b(0,0))>(std::forward<TypeA>(a).base + std::forward<TypeB>(b).base);
    }

    // matrix + scalar
    template<typename TypeA, typename TypeB, typename enable = typename std::enable_if<
        meta::is_matrix<TypeA>::value && meta::is_scalar<TypeB>::value
    >::type>
    auto operator - (TypeA&& a, const TypeB& b) -> mat<decltype(a(0,0)+b)> {
        return mat<decltype(a(0,0)+b)>(std::forward<TypeA>(a).base + b);
    }

    // scalar + matrix
    template<typename TypeA, typename TypeB, typename enable = typename std::enable_if<
        meta::is_matrix<TypeB>::value && meta::is_scalar<TypeA>::value
    >::type>
    auto operator - (const TypeA& a, TypeB&& b) -> mat<decltype(a+b(0,0))> {
        return mat<decltype(a+b(0,0))>(a + std::forward<TypeB>(b).base);
    }

    // matrix - matrix
    template<typename TypeA, typename TypeB, typename enable = typename std::enable_if<
        meta::is_matrix<TypeA>::value && meta::is_matrix<TypeB>::value
    >::type>
    auto operator - (TypeA&& a, TypeB&& b) -> mat<decltype(a(0,0)-b(0,0))> {
        return mat<decltype(a(0,0)-b(0,0))>(std::forward<TypeA>(a).base - std::forward<TypeB>(b).base);
    }

    // matrix - scalar
    template<typename TypeA, typename TypeB, typename enable = typename std::enable_if<
        meta::is_matrix<TypeA>::value && meta::is_scalar<TypeB>::value
    >::type>
    auto operator - (TypeA&& a, const TypeB& b) -> mat<decltype(a(0,0)-b)> {
        return mat<decltype(a(0,0)-b)>(std::forward<TypeA>(a).base - b);
    }

    // scalar - matrix
    template<typename TypeA, typename TypeB, typename enable = typename std::enable_if<
        meta::is_matrix<TypeB>::value && meta::is_scalar<TypeA>::value
    >::type>
    auto operator - (const TypeA& a, TypeB&& b) -> mat<decltype(a-b(0,0))> {
        return mat<decltype(a-b(0,0))>(a - std::forward<TypeB>(b).base);
    }
}
}
