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
        vif_check(v.dims[0] == v.dims[1], "cannot get diagnoal of non-square matrix (got ",
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
        vif_check(v.dims[0] == v.dims[1], "cannot get diagnoal of non-square matrix (got ",
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

    // Source of this class is adapted from:
    // http://www.mymathlib.com/matrices/linearsystems/crout.html
    struct decompose_lu {
        // Options
        bool no_pivot = false;

        // Outputs
        mat2d lu;
        vec1u ipiv;
        bool bad = false;
        uint_t ns = 0;

    public:
        uint_t size() const {
            return ipiv.size();
        }

        template<typename T, typename enable = typename std::enable_if<
            meta::is_matrix<T>::value
        >::type>
        bool decompose(T alpha) {
            // LU decomposition

            vif_check(alpha.dims[0] == alpha.dims[1], "cannot do LU decomposition of a non "
                "square matrix (", alpha.dims, ")");

            const uint_t n = alpha.dims[0];

            bad = false;
            lu = std::move(alpha);
            ipiv = indgen(n);
            ns = 0;

            for (uint_t k : range(n)) {
                // Find pivot
                double akk = abs(lu.safe(k,k));
                if (!no_pivot) {
                    for (uint_t j : range(k+1, n)) {
                        double ajk = abs(lu.safe(j,k));
                        if (akk < ajk) {
                            akk = ajk;
                            ipiv.safe[k] = j;
                        }
                    }

                    // Apply pivot
                    if (ipiv.safe[k] != k) {
                        ++ns;

                        for (uint_t j : range(n)) {
                            std::swap(lu.safe(ipiv.safe[k],j), lu.safe(k,j));
                        }
                    }
                }

                // If the matrix is singular, return an error
                if (lu.safe(k,k) == 0.0) {
                    bad = true;
                }

                // Find upper elements of row
                for (uint_t j : range(k+1, n)) {
                    lu.safe(k,j) /= lu.safe(k,k);
                }

                // Update remaining matrix elements
                for (uint_t i : range(k+1, n))
                for (uint_t j : range(k+1, n)) {
                    lu.safe(i,j) -= lu.safe(i,k)*lu.safe(k,j);
                }
            }

            return !bad;
        }

        vec1d solve(vec1d x) const {
            vif_check(lu.dims[0] == x.dims[0], "matrix and vector must have the same "
                "dimensions (got ", lu.dims[0], " and ", x.dims[0], ")");

            vec1d r(x.size());
            const uint_t n = x.size();

            // Solve L*y = x
            for (uint_t k : range(n)) {
                if (!no_pivot) {
                    // Apply pivot
                    uint_t pk = ipiv.safe[k];
                    if (pk != k) {
                        std::swap(x.safe[k], x.safe[pk]);
                    }
                }

                r.safe[k] = x.safe[k];
                for (uint_t i : range(k)) {
                    r.safe[k] -= r.safe[i]*lu.safe(k,i);
                }

                r.safe[k] /= lu.safe(k,k);
            }

            // Solve U*z = y
            uint_t k = n;
            while (k > 0) {
                --k;

                for (uint_t i : range(k+1, n)) {
                    r.safe[k] -= r.safe[i]*lu.safe(k,i);
                }
            }

            return r;
        }

        mat2d invert() const {
            mat2d inv(lu.dims);
            const uint_t n = lu.dims[0];

            vec1d r(n);
            for (uint_t c : range(n)) {
                uint_t i1 = c;

                // Solve L*y = x
                for (uint_t k : range(n)) {
                    if (!no_pivot) {
                        // Apply pivot
                        uint_t pk = ipiv.safe[k];
                        if (k == i1) {
                            i1 = pk;
                        } else if (pk == i1) {
                            i1 = k;
                        }

                        r.safe[k] = (i1 == k);
                    } else {
                        r.safe[k] = 1.0;
                    }

                    for (uint_t i : range(k)) {
                        r.safe[k] -= r.safe[i]*lu.safe(k,i);
                    }

                    r.safe[k] /= lu.safe(k,k);
                }

                // Solve U*z = y
                uint_t k = n;
                while (k > 0) {
                    --k;

                    for (uint_t i : range(k+1, n)) {
                        r.safe[k] -= r.safe[i]*lu.safe(k,i);
                    }

                    inv.safe(k,c) = r.safe[k];
                }
            }

            return inv;
        }

        double determinant() const {
            double d = (ns % 2 == 0 ? 1.0 : -1.0);
            const uint_t n = lu.dims[0];
            for (uint_t i : range(n)) {
                d *= lu.safe(i,i);
            }

            return d;
        }
    };

    // Source of this class is adapted from:
    // https://rosettacode.org/wiki/Cholesky_decomposition#C
    struct decompose_cholesky {
        // Outputs
        mat2d l;
        bool bad = false;

    public:
        uint_t size() const {
            return l.dims[0];
        }

        template<typename T, typename enable = typename std::enable_if<
            meta::is_matrix<T>::value
        >::type>
        bool decompose(T alpha) {
            // Cholesky decomposition
            vif_check(alpha.dims[0] == alpha.dims[1], "cannot do Cholesky decomposition of a non "
                "square matrix (", alpha.dims, ")");

            bad = false;
            l = std::move(alpha);

        #ifdef NO_LAPACK
            const uint_t n = alpha.dims[0];
            for (uint_t i : range(n)) {
                for (uint_t j : range(i+1)) {
                    double s = 0.0;
                    for (uint_t k : range(j)) {
                        s += l.safe(i,k)*l.safe(j,k);
                    }

                    if (i == j) {
                        s = l.safe(i,i) - s;
                        if (s <= 0.0) {
                            bad = true;
                            return false;
                        }
                        l.safe(i,i) = sqrt(s);
                    } else {
                        l.safe(i,j) = (l.safe(i,j) - s)/l.safe(j,j);
                    }
                }

                for (uint_t j : range(i+1, n)) {
                    l.safe(i,j) = 0.0;
                }
            }
        #else
            char uplo = 'U';
            int n = l.dims[0];
            int lda = n;
            int info;

            lapack::dpotrf_(&uplo, &n, l.raw_data(), &lda, &info);

            for (uint_t i : range(n))
            for (uint_t j : range(i+1, n)) {
                l.safe(i,j) = 0.0;
            }

            bad = info != 0;
        #endif

            return !bad;
        }

        void substitute_forward_inplace(vec1d& x) const {
            vif_check(l.dims[0] == x.dims[0], "matrix and vector must have the same "
                "dimensions (got ", l.dims[0], " and ", x.dims[0], ")");

            const uint_t n = x.size();

            // Solve L*y = x
            for (uint_t k : range(n)) {
                for (uint_t i : range(k)) {
                    x.safe[k] -= x.safe[i]*l.safe(k,i);
                }

                x.safe[k] /= l.safe(k,k);
            }
        }

        vec1d substitute_forward(vec1d x) const {
            substitute_forward_inplace(x);
            return x;
        }

        void substitute_backward_inplace(vec1d& x) const {
            vif_check(l.dims[0] == x.dims[0], "matrix and vector must have the same "
                "dimensions (got ", l.dims[0], " and ", x.dims[0], ")");

            const uint_t n = x.size();

            // Solve L^t*y = x
            uint_t k = n;
            while (k > 0) {
                --k;

                for (uint_t i : range(k+1, n)) {
                    x.safe[k] -= x.safe[i]*l.safe(i,k);
                }

                x.safe[k] /= l.safe(k,k);
            }
        }

        vec1d substitute_backward(vec1d x) const {
            substitute_backward_inplace(x);
            return x;
        }

        void solve_inplace(vec1d& x) const {
            vif_check(l.dims[0] == x.dims[0], "matrix and vector must have the same "
                "dimensions (got ", l.dims[0], " and ", x.dims[0], ")");

        #ifdef NO_LAPACK
            substitute_forward_inplace(x);
            substitute_backward_inplace(x);
        #else
            char uplo = 'U';
            int n = l.dims[0];
            int nrhs = 1;
            int lda = n;
            int ldb = n;
            int info;

            lapack::dpotrs_(&uplo, &n, &nrhs, l.raw_data(), &lda, x.raw_data(), &ldb, &info);
        #endif
        }

        vec1d solve(vec1d x) const {
            solve_inplace(x);
            return x;
        }

        mat2d invert() const {
        #ifdef NO_LAPACK
            mat2d inv(l.dims);
            const uint_t n = l.dims[0];

            vec1d x(n);
            for (uint_t c : range(n)) {
                x.safe[c] = 1.0;
                substitute_forward_inplace(x);
                substitute_backward_inplace(x);
                for (uint_t i : range(x)) {
                    inv.safe(c,i) = x.safe[i];
                    x.safe[i] = 0.0;
                }
            }
        #else
            mat2d inv = l;
            char uplo = 'U';
            int n = l.dims[0];
            int lda = n;
            int info;

            lapack::dpotri_(&uplo, &n, inv.raw_data(), &lda, &info);

            for (uint_t i : range(n))
            for (uint_t j : range(i+1, n)) {
                inv.safe(i,j) = inv.safe(j,i);
            }
        #endif

            return inv;
        }

        mat2d lower_inverse() const {
        #ifdef NO_LAPACK
            mat2d inv(l.dims);
            const uint_t n = l.dims[0];

            vec1d x(n);
            for (uint_t c : range(n)) {
                x.safe[c] = 1.0;
                substitute_backward_inplace(x);
                for (uint_t i : range(x)) {
                    inv.safe(c,i) = x.safe[i];
                    x.safe[i] = 0.0;
                }
            }
        #else
            mat2d inv = l;
            char uplo = 'U';
            char diag = 'N';
            int n = l.dims[0];
            int lda = n;
            int info = 0;

            lapack::dtrtri_(&uplo, &diag, &n, inv.raw_data(), &lda, &info);
        #endif

            return inv;
        }

        double determinant() const {
            double d = 1.0;
            const uint_t n = l.dims[0];
            for (uint_t i : range(n)) {
                d *= l.safe(i,i);
            }

            return d;
        }

        double log_determinant() const {
            double d = 0.0;
            const uint_t n = l.dims[0];
            for (uint_t i : range(n)) {
                d += log(l.safe(i,i));
            }

            return d;
        }
    };


    template<typename Type, typename enable = typename std::enable_if<
        meta::is_matrix<Type>::value
    >::type>
    auto determinant(const Type& v) -> decltype(+v(0,0)) {
        vif_check(v.dims[0] == v.dims[1], "cannot compute determinant of non-square matrix (got ",
            v.dims, ")");

        decompose_lu d;
        if (!d.decompose(v)) {
            return 0.0;
        }

        return d.determinant();
    }

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
    bool solve(const T& alpha, const vec1d& beta, vec1d& res);

    template<typename T, typename enable = typename std::enable_if<
        meta::is_matrix<T>::value && std::is_same<meta::data_type_t<T>, double>::value &&
        !meta::is_view<T>::value
    >::type>
    bool inplace_solve(T& alpha, vec1d& beta) {
    #ifdef NO_LAPACK
        vec1d cbeta = beta;
        return solve(alpha, cbeta, beta);
    #else
        vif_check(alpha.dims[0] == alpha.dims[1], "cannot invert a non square matrix (",
            alpha.dims, ")");
        vif_check(alpha.dims[0] == beta.dims[0], "matrix and vector must have the same dimensions (",
            "got ", alpha.dims[0], " and ", beta.dims[0], ")");

        int n = alpha.dims[0];
        int nrhs = 1;
        int lda = n, ldb = n;
        int info;
        vec<1,int> ipiv(n);

        lapack::dgesv_(&n, &nrhs, alpha.raw_data(), &lda, ipiv.raw_data(), beta.raw_data(),
            &ldb, &info);
        if (info != 0) {
            return false;
        }

        return true;
    #endif
    }

    template<typename T, typename enable>
    bool solve(const T& alpha, const vec1d& beta, vec1d& res) {
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
        return inplace_solve(a, res);
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
