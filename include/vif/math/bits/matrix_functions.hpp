#ifndef VIF_INCLUDING_MATH_MATRIX_BITS
#error this file is not meant to be included separately, include "vif/math/matrix.hpp" instead
#endif

namespace vif {
namespace matrix {
    template<typename Type, typename enable = typename std::enable_if<
        meta::is_matrix<Type>::value && !std::is_reference<Type>::value
    >::type>
    mat<typename meta::data_type<typename std::decay<Type>::type>::type> transpose(const Type& v) {
        using rtype = mat<typename meta::data_type<typename std::decay<Type>::type>::type>;
        return rtype(vif::transpose(v.base));
    }

    template<typename Type, typename enable = typename std::enable_if<
        meta::is_matrix<Type>::value && !std::is_reference<Type>::value
    >::type>
    auto diagonal(Type&& v) -> decltype(v(_,0).concretise()) {
        vif_check(v.dims[0] == v.dims[1], "can only be called on square matrix (got ",
            v.dims, ")");

        decltype(v(_,0).concretise()) d(v.dims[0]);
        for (uint_t i : range(d)) {
            d.safe[i] = v.safe(i,i);
        }

        return d;
    }

    template<typename Type, typename enable = typename std::enable_if<
        meta::is_matrix<Type>::value
    >::type>
    auto diagonal(Type& v) -> decltype(v(_,0)) {
        vif_check(v.dims[0] == v.dims[1], "can only be called on square matrix (got ",
            v.dims, ")");

        decltype(v(_,0)) d(impl::vec_ref_tag, impl::vec_access::get_parent(v.base));
        d.dims[0] = v.dims[0];
        d.resize();
        for (uint_t i : range(d)) {
            d.data[i] = impl::ptr<typename meta::data_type<Type>::type>(v.safe(i,i));
        }

        return d;
    }

    template<typename Type = double>
    mat<Type> make_identity(uint_t dim) {
        mat<Type> m(dim, dim);
        diagonal(m) = 1;
        return m;
    }

    template<typename TX, typename TY>
    auto make_scale(const TX& sx, const TY& sy) -> mat<decltype(sx*sy)> {
        mat<decltype(sx*sy)> m(3, 3);
        m.safe(0,0) = sx;
        m.safe(1,1) = sy;
        m.safe(2,2) = 1;
        return m;
    }

    template<typename T>
    mat<T> make_scale(const T& s) {
        mat<T> m(3, 3);
        m.safe(0,0) = s;
        m.safe(1,1) = s;
        m.safe(2,2) = 1;
        return m;
    }

    template<typename TX, typename TY>
    auto make_translation(const TX& tx, const TY& ty) -> mat<decltype(tx*ty)> {
        mat<decltype(tx*ty)> m(3, 3);
        diagonal(m) = 1;
        m.safe(0,2) = tx;
        m.safe(1,2) = ty;
        return m;
    }

    template<typename A>
    auto make_rotation(const A& a) -> mat<decltype(cos(a))> {
        mat<decltype(cos(a))> m(3, 3);
        auto ca = cos(a), sa = sin(a);
        m.safe(0,0) = m.safe(1,1) = ca;
        m.safe(0,1) = -sa;
        m.safe(1,0) = sa;
        m.safe(2,2) = 1;
        return m;
    }

    template<typename Type, typename enable = typename std::enable_if<
        meta::is_matrix<Type>::value
    >::type>
    void inplace_symmetrize(Type& alpha) {
        vif_check(alpha.dims[0] == alpha.dims[1], "cannot symmetrize a non square matrix (",
            alpha.dims, ")");

        for (uint_t i : range(alpha.dims[0]))
        for (uint_t j : range(i+1, alpha.dims[0])) {
            alpha.safe(i,j) = alpha.safe(j,i);
        }
    }

    struct decompose_lu {
        mat2d lu;
        vec1u ipiv;
        mutable vec1d y;
        bool bad = false;

    public:
        template<typename T, typename enable = typename std::enable_if<
            meta::is_matrix<T>::value
        >::type>
        bool decompose(const T& alpha) {
            // LU decomposition

            vif_check(alpha.dims[0] == alpha.dims[1], "cannot do LU decomposition of a non "
                "square matrix (", alpha.dims, ")");

            const uint_t n = alpha.dims[0];

            y.resize(n);
            lu.resize(alpha.dims);
            ipiv = indgen<uint_t>(n);

            // Find pivot
            for (uint_t i : range(n)) {
                double aii = abs(alpha.safe(i,i));
                for (uint_t k : range(i+1, n)) {
                    double aki = abs(alpha.safe(k,i));
                    if (aii < aki) {
                        aii = aki;
                        std::swap(ipiv.safe[i], ipiv.safe[k]);
                    }
                }
            }

            // Find L and U
            for (uint_t k : range(n)) {
                for (uint_t i : range(k, n)) {
                    double sum = 0.0;
                    for (uint_t p : range(k)) {
                        sum += lu.safe(i,p)*lu.safe(p,k);
                    }

                    lu.safe(i,k) = alpha.safe(ipiv.safe[i],k) - sum;
                }

                if (lu.safe(k,k) == 0.0) {
                    bad = true;
                }

                for (uint_t i : range(k+1, n)) {
                    double sum = 0.0;
                    for (uint_t p : range(k)) {
                        sum += lu.safe(k,p)*lu.safe(p,i);
                    }

                    lu.safe(k,i) = (alpha.safe(ipiv.safe[k],i) - sum)/lu.safe(k,k);
                }
            }

            return !bad;
        }

        mat2d invert() const {
            mat2d inv(lu.dims);
            const uint_t n = lu.dims[0];

            for (uint_t c : range(n)) {
                // Solve L*y = x
                for (uint_t i : range(n)) {
                    double sum = 0.0;
                    for (uint_t p : range(i)) {
                        sum += lu.safe(i,p)*y.safe[p];
                    }

                    y.safe[i] = ((c == ipiv.safe[i] ? 1.0 : 0.0) - sum)/lu.safe(i,i);
                }

                // Solve U*z = y
                {
                    uint_t i = n;
                    while (i > 0) {
                        --i;

                        double sum = 0.0;
                        uint_t p = n;
                        while (p > i+1) {
                            --p;
                            sum += lu.safe(i,p)*inv.safe(p,c);
                        }

                        inv.safe(i,c) = y.safe[i] - sum;
                    }
                }
            }

            return inv;
        }

        vec1d solve(const vec1d& x) const {
            vif_check(lu.dims[0] == x.dims[0], "matrix and vector must have the same "
                "dimensions (got ", lu.dims[0], " and ", x.dims[0], ")");

            vec1d r(x.size());
            const uint_t n = x.size();

            // Solve L*y = x
            for (uint_t i : range(n)) {
                double sum = 0.0;
                for (uint_t p : range(i)) {
                    sum += lu.safe(i,p)*y.safe[p];
                }

                y.safe[i] = (x.safe[ipiv.safe[i]] - sum)/lu.safe(i,i);
            }

            // Solve U*z = y
            {
                uint_t i = n;
                while (i > 0) {
                    --i;

                    double sum = 0.0;
                    uint_t p = n;
                    while (p > i+1) {
                        --p;
                        sum += lu.safe(i,p)*r.safe[p];
                    }

                    r.safe[i] = y.safe[i] - sum;
                }
            }

            return r;
        }
    };

    template<typename T, typename enable = typename std::enable_if<
        meta::is_matrix<T>::value && std::is_same<meta::data_type_t<T>, double>::value &&
        !meta::is_view<T>::value
    >::type>
    bool invert(const T& a, T& i);

    template<typename T, typename enable = typename std::enable_if<
        meta::is_matrix<T>::value && std::is_same<meta::data_type_t<T>, double>::value &&
        !meta::is_view<T>::value
    >::type>
    bool inplace_invert(T& i) {
    #ifdef NO_LAPACK
        mat2d a = i;
        return invert(a, i);
    #else
        vif_check(i.dims[0] == i.dims[1], "cannot invert a non square matrix (", i.dims, ")");

        int n = i.dims[0];
        int lda = n;
        int info;

        vec<1,int> ipiv(n);
        lapack::dgetrf_(&n, &n, i.raw_data(), &lda, ipiv.raw_data(), &info);
        if (info < 0) {
            return false;
        }

        vec1d work(n);
        int lw = n;
        lapack::dgetri_(&n, i.raw_data(), &lda, ipiv.raw_data(), work.raw_data(), &lw, &info);
        if (info != 0) {
            return false;
        }

        return true;
    #endif
    }

    template<typename T, typename enable>
    bool invert(const T& a, T& i) {
    #ifdef NO_LAPACK
        decompose_lu d;
        if (!d.decompose(a)) {
            return false;
        }

        i = d.invert();
        return true;
    #else
        i = a;
        return inplace_invert(i);
    #endif
    }

    template<typename T, typename enable = typename std::enable_if<
        meta::is_matrix<T>::value && std::is_same<meta::data_type_t<T>, double>::value &&
        !meta::is_view<T>::value
    >::type>
    bool invert_symmetric(const T& a, T& i);

    template<typename T, typename enable = typename std::enable_if<
        meta::is_matrix<T>::value && std::is_same<meta::data_type_t<T>, double>::value &&
        !meta::is_view<T>::value
    >::type>
    bool inplace_invert_symmetric(T& i) {
    #ifdef NO_LAPACK
        mat2d a = i;
        return invert_symmetric(a, i);
    #else
        vif_check(i.dims[0] == i.dims[1], "cannot invert a non square matrix (", i.dims, ")");

        char uplo = 'U';
        int n = i.dims[0];
        int lda = n;
        int info;

        int lw = n*64;
        // Note: the optimal value for lw is n*nb, where nb is the optimal block size
        // This value can be obtained using ilaenv_, but 64 should be plenty enough, according to
        // the Lapack User Guide.

        vec1d work(lw);
        vec<1,int> ipiv(n);

        lapack::dsytrf_(&uplo, &n, i.raw_data(), &lda, ipiv.raw_data(), work.raw_data(), &lw, &info);
        if (info < 0) {
            return false;
        }

        lapack::dsytri_(&uplo, &n, i.raw_data(), &lda, ipiv.raw_data(), work.raw_data(), &info);
        if (info != 0) {
            return false;
        }

        return true;
    #endif
    }

    template<typename T, typename enable>
    bool invert_symmetric(const T& a, T& i) {
    #ifdef NO_LAPACK
        decompose_lu d;
        if (!d.decompose(a)) {
            return false;
        }

        i = d.invert();
        return true;
    #else
        i = a;
        return inplace_invert_symmetric(i);
    #endif
    }

    template<typename T, typename enable = typename std::enable_if<
        meta::is_matrix<T>::value && std::is_same<meta::data_type_t<T>, double>::value &&
        !meta::is_view<T>::value
    >::type>
    bool solve_symmetric(const T& alpha, const vec1d& beta, vec1d& res);

    template<typename T, typename enable = typename std::enable_if<
        meta::is_matrix<T>::value && std::is_same<meta::data_type_t<T>, double>::value &&
        !meta::is_view<T>::value
    >::type>
    bool inplace_solve_symmetric(T& alpha, vec1d& beta) {
    #ifdef NO_LAPACK
        vec1d cbeta = beta;
        return solve_symmetric(alpha, cbeta, beta);
    #else
        vif_check(alpha.dims[0] == alpha.dims[1], "cannot invert a non square matrix (",
            alpha.dims, ")");
        vif_check(alpha.dims[0] == beta.dims[0], "matrix and vector must have the same dimensions (",
            "got ", alpha.dims[0], " and ", beta.dims[0], ")");

        char uplo = 'U';
        int n = alpha.dims[0];
        int nrhs = 1;
        int lda = n, ldb = n;
        int info;

        int lw = n*64;
        // Note: the optimal value for lw is n*nb, where nb is the optimal block size
        // This value can be obtained using ilaenv_, but 64 should be plenty enough, according to
        // the Lapack User Guide.

        vec1d work(lw);
        vec<1,int> ipiv(n);

        lapack::dsysv_(&uplo, &n, &nrhs, alpha.raw_data(), &lda, ipiv.raw_data(), beta.raw_data(),
            &ldb, work.raw_data(), &lw, &info);
        if (info != 0) {
            return false;
        }

        return true;
    #endif
    }

    template<typename T, typename enable>
    bool solve_symmetric(const T& alpha, const vec1d& beta, vec1d& res) {
    #ifdef NO_LAPACK
        decompose_lu d;
        if (!d.decompose(alpha)) {
            return false;
        }

        res = d.solve(beta);
        return true;
    #else
        mat2d a = alpha;
        res = beta;
        return inplace_solve_symmetric(a, res);
    #endif
    }

    template<typename T, typename enable = typename std::enable_if<
        meta::is_matrix<T>::value && std::is_same<meta::data_type_t<T>, double>::value &&
        !meta::is_view<T>::value
    >::type>
    bool inplace_eigen_symmetric(T& a, vec1d& vals) {
    #ifdef NO_LAPACK
        static_assert(!std::is_same<T,T>::value, "LAPACK support has been disabled, "
            "please enable LAPACK to use this function");
        return false;
    #else
        vif_check(a.dims[0] == a.dims[1], "cannot invert a non square matrix (",
            a.dims, ")");

        char jobz = 'V';
        char uplo = 'U';
        int n = a.dims[0];
        int lda = n;
        int info;

        vals.resize(n);

        int lw = n*64;
        // Note: the optimal value for lw is n*nb, where nb is the optimal block size
        // This value can be obtained using ilaenv_, but 64 should be plenty enough, according to
        // the Lapack User Guide.

        vec1d work(lw);

        lapack::dsyev_(&jobz, &uplo, &n, a.raw_data(), &lda, vals.raw_data(), work.raw_data(),
            &lw, &info);
        if (info != 0) {
            return false;
        }

        // Eigen vectors are now stored in 'a' with the following layout:
        // v0 = a(0,_), v1 = a(1,_), ...
        // each corresponding to the eigen values given in 'vals'

        return true;
    #endif
    }

    template<typename T, typename enable = typename std::enable_if<
        meta::is_matrix<T>::value && std::is_same<meta::data_type_t<T>, double>::value &&
        !meta::is_view<T>::value
    >::type>
    bool eigen_symmetric(const T& a, vec1d& vals, T& vecs) {
        vecs = a;
        return inplace_eigen_symmetric(vecs, vals);
    }
}
}
