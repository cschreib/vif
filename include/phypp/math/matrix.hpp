#ifndef PHYPP_MATH_MATRIX_HPP
#define PHYPP_MATH_MATRIX_HPP

#ifndef NO_LAPACK
#include "phypp/math/lapack.hpp"
#endif
#include "phypp/core/vec.hpp"
#include "phypp/core/error.hpp"
#include "phypp/core/range.hpp"
#include "phypp/math/base.hpp"

namespace phypp {
namespace matrix {
    using phypp::transpose;

    template<typename TypeA, typename TypeB>
    auto product(const vec<2,TypeA>& a, const vec<2,TypeB>& b) -> vec<2,decltype(a(0,0)*b(0,0))> {
        phypp_check(a.dims[1] == b.dims[0], "incompatible dimensions in matrix-matrix multiplication "
            "(", a.dims, " x ", b.dims, ")");

        const uint_t o = a.dims[1];

        using ntype_t = decltype(a(0,0)*b(0,0));
        const uint_t n = a.dims[0];
        const uint_t m = b.dims[1];

        vec<2,ntype_t> r(n,m);
        for (uint_t i : range(n))
        for (uint_t j : range(m))
        for (uint_t k : range(o)) {
            r.safe(i,j) += a.safe(i,k)*b.safe(k,j);
        }

        return r;
    }

    template<typename TypeA, typename TypeB>
    auto product(const vec<2,TypeA>& a, const vec<1,TypeB>& b) -> vec<1,decltype(a(0,0)*b(0,0))> {
        phypp_check(a.dims[1] == b.dims[0], "incompatible dimensions in matrix-vector multiplication "
            "(", a.dims, " x ", b.dims, ")");

        const uint_t o = a.dims[1];

        using ntype_t = decltype(a(0,0)*b(0,0));
        const uint_t n = a.dims[0];

        vec<1,ntype_t> r(n);
        for (uint_t i : range(n))
        for (uint_t k : range(o)) {
            r.safe(i) += a.safe(i,k)*b.safe(k);
        }

        return r;
    }

    template<typename TypeA, typename TypeB>
    auto product(const vec<1,TypeB>& b, const vec<2,TypeA>& a) -> vec<1,decltype(a(0,0)*b(0,0))> {
        phypp_check(a.dims[1] == b.dims[0], "incompatible dimensions in vector-matrix multiplication "
            "(", a.dims, " x ", b.dims, ")");

        const uint_t o = a.dims[0];

        using ntype_t = decltype(a(0,0)*b(0,0));
        const uint_t n = a.dims[1];

        vec<1,ntype_t> r = arr<ntype_t>(n);
        for (uint_t i : range(n))
        for (uint_t k : range(o)) {
            r.safe(i) += b.safe(k)*a.safe(k,i);
        }

        return r;
    }

    template<typename Type>
    auto diagonal(vec<2,Type>&& v) -> decltype(v(_,0).concretise()) {
        phypp_check(v.dims[0] == v.dims[1], "can only be called on square matrix (got ",
            v.dims, ")");

        decltype(v(_,0).concretise()) d(v.dims[0]);
        for (uint_t i : range(d)) {
            d.safe[i] = v.safe(i,i);
        }

        return d;
    }

    template<typename Type>
    auto diagonal(const vec<2,Type>& v) -> decltype(v(_,0)) {
        phypp_check(v.dims[0] == v.dims[1], "can only be called on square matrix (got ",
            v.dims, ")");

        decltype(v(_,0)) d(impl::vec_ref_tag, impl::vec_access::get_parent(v));
        d.dims[0] = v.dims[0];
        d.resize();
        for (uint_t i : range(d)) {
            d.safe[i] = impl::ptr<Type>(v.safe(i,i));
        }

        return d;
    }

    template<typename Type>
    auto diagonal(vec<2,Type>& v) -> decltype(v(_,0)) {
        phypp_check(v.dims[0] == v.dims[1], "can only be called on square matrix (got ",
            v.dims, ")");

        decltype(v(_,0)) d(impl::vec_ref_tag, impl::vec_access::get_parent(v));
        d.dims[0] = v.dims[0];
        d.resize();
        for (uint_t i : range(d)) {
            d.data[i] = impl::ptr<Type>(v.safe(i,i));
        }

        return d;
    }

    template<typename Type = double>
    vec<2,Type> make_identity(uint_t dim) {
        vec<2,Type> m(dim, dim);
        diagonal(m) = 1;
        return m;
    }

    template<typename TX, typename TY>
    auto make_scale(const TX& sx, const TY& sy) -> vec<2,decltype(sx*sy)> {
        vec<2,decltype(sx*sy)> m(3, 3);
        m.safe(0,0) = sx;
        m.safe(1,1) = sy;
        m.safe(2,2) = 1;
        return m;
    }

    template<typename T>
    vec<2,T> make_scale(const T& s) {
        vec<2,T> m(3, 3);
        m.safe(0,0) = s;
        m.safe(1,1) = s;
        m.safe(2,2) = 1;
        return m;
    }

    template<typename TX, typename TY>
    auto make_translation(const TX& tx, const TY& ty) -> vec<2,decltype(tx*ty)> {
        vec<2,decltype(tx*ty)> m(3, 3);
        diagonal(m) = 1;
        m.safe(0,2) = tx;
        m.safe(1,2) = ty;
        return m;
    }

    template<typename A>
    auto make_rotation(const A& a) -> vec<2,decltype(cos(a))> {
        vec<2,decltype(cos(a))> m(3, 3);
        auto ca = cos(a), sa = sin(a);
        m.safe(0,0) = m.safe(1,1) = ca;
        m.safe(0,1) = -sa;
        m.safe(1,0) = sa;
        m.safe(2,2) = 1;
        return m;
    }

    template<typename Type>
    void symmetrize(vec<2,Type>& alpha) {
        phypp_check(alpha.dims[0] == alpha.dims[1], "cannot symmetrize a non square matrix (",
            alpha.dims, ")");

        for (uint_t i : range(alpha.dims[0]))
        for (uint_t j : range(i+1, alpha.dims[0])) {
            alpha.safe(i,j) = alpha.safe(j,i);
        }
    }

    template<typename TX, typename TY>
    auto make_point(const TX& x, const TY& y) -> vec<1,decltype(x*y)> {
        return vec<1,decltype(x*y)>{x, y, 1};
    }

    template<typename Dummy = void>
    bool inplace_invert(vec2d& i) {
    #ifdef NO_LAPACK
        static_assert(!std::is_same<Dummy,Dummy>::value, "LAPACK support has been disabled, "
            "please enable LAPACK to use this function");
        return false;
    #else
        phypp_check(i.dims[0] == i.dims[1], "cannot invert a non square matrix (", i.dims, ")");

        int n = i.dims[0];
        int lda = n;
        int info;

        vec<1,int> ipiv(n);
        lapack::dgetrf_(&n, &n, i.data.data(), &lda, ipiv.data.data(), &info);
        if (info < 0) {
            return false;
        }

        vec1d work(n);
        int lw = n;
        lapack::dgetri_(&n, i.data.data(), &lda, ipiv.data.data(), work.data.data(), &lw, &info);
        if (info != 0) {
            return false;
        }

        return true;
    #endif
    }

    template<typename TypeA>
    bool invert(const vec<2,TypeA>& a, vec2d& i) {
        i = a;
        return inplace_invert<TypeA>(i);
    }

    template<typename Dummy = void>
    bool inplace_invert_symmetric(vec2d& i) {
    #ifdef NO_LAPACK
        static_assert(!std::is_same<Dummy,Dummy>::value, "LAPACK support has been disabled, "
            "please enable LAPACK to use this function");
        return false;
    #else
        phypp_check(i.dims[0] == i.dims[1], "cannot invert a non square matrix (", i.dims, ")");

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

        lapack::dsytrf_(&uplo, &n, i.data.data(), &lda, ipiv.data.data(), work.data.data(), &lw, &info);
        if (info < 0) {
            return false;
        }

        lapack::dsytri_(&uplo, &n, i.data.data(), &lda, ipiv.data.data(), work.data.data(), &info);
        if (info != 0) {
            return false;
        }

        return true;
    #endif
    }

    template<typename TypeA>
    bool invert_symmetric(const vec<2,TypeA>& a, vec2d& i) {
        i = a;
        return inplace_invert_symmetric<TypeA>(i);
    }

    template<typename Dummy = void>
    bool inplace_solve_symmetric(vec2d& alpha, vec1d& beta) {
    #ifdef NO_LAPACK
        static_assert(!std::is_same<Dummy,Dummy>::value, "LAPACK support has been disabled, "
            "please enable LAPACK to use this function");
        return false;
    #else
        phypp_check(alpha.dims[0] == alpha.dims[1], "cannot invert a non square matrix (",
            alpha.dims, ")");
        phypp_check(alpha.dims[0] == beta.dims[0], "matrix and vector must have the same dimensions (",
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

        lapack::dsysv_(&uplo, &n, &nrhs, alpha.data.data(), &lda, ipiv.data.data(), beta.data.data(),
            &ldb, work.data.data(), &lw, &info);
        if (info != 0) {
            return false;
        }

        return true;
    #endif
    }

    template<typename Dummy = void>
    bool solve_symmetric(const vec2d& alpha, const vec1d& beta, vec1d& res) {
        vec2d a = alpha;
        res = beta;
        return inplace_solve_symmetric<Dummy>(a, res);
    }

    template<typename Dummy = void>
    bool inplace_eigen_symmetric(vec2d& a, vec1d& vals) {
    #ifdef NO_LAPACK
        static_assert(!std::is_same<Dummy,Dummy>::value, "LAPACK support has been disabled, "
            "please enable LAPACK to use this function");
        return false;
    #else
        phypp_check(a.dims[0] == a.dims[1], "cannot invert a non square matrix (",
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

        lapack::dsyev_(&jobz, &uplo, &n, a.data.data(), &lda, vals.data.data(), work.data.data(),
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

    template<typename Dummy = void>
    bool eigen_symmetric(const vec2d& a, vec1d& vals, vec2d& vecs) {
        vecs = a;
        return inplace_eigen_symmetric<Dummy>(vecs, vals);
    }

    struct linear_solver {
        vec2d alpha;
        vec2d l, u;
        vec1u ipiv;
        mutable vec1d y, z;

        void decompose() {
            const uint_t n = alpha.dims[0];

            y.resize(n);
            z.resize(n);
            l.resize(alpha.dims);
            u.resize(alpha.dims);
            ipiv = uindgen(n);

            // Find and apply pivot
            for (uint_t i : range(n)) {
                double aii = abs(alpha.safe(i,i));
                for (uint_t k : range(i+1, n)) {
                    double aki = abs(alpha.safe(k,i));
                    if (aii < aki) {
                        aii = aki;
                        std::swap(ipiv.safe[i], ipiv.safe[k]);
                        for (uint_t j : range(n)) {
                            std::swap(alpha.safe(i,j), alpha.safe(k,j));
                        }
                    }
                }
            }

            // Find L and U
            for (uint_t k : range(n)) {
                u.safe(k,k) = 1.0;
                for (uint_t i : range(k, n)) {
                    double sum = 0.0;
                    for (uint_t p : range(k)) {
                        sum += l.safe(i,p)*u.safe(p,k);
                    }

                    l.safe(i,k) = alpha.safe(i,k) - sum;
                }

                for (uint_t i : range(k+1, n)) {
                    double sum = 0.0;
                    for (uint_t p : range(k)) {
                        sum += l.safe(k,p)*u.safe(p,i);
                    }

                    u.safe(k,i) = (alpha.safe(k,i) - sum)/l.safe(k,k);
                }
            }
        }

        vec1d solve(const vec1d& x) const {
            vec1d r(x.size());

            const uint_t n = x.size();

            // Solve L*y = x
            for (uint_t i : range(n)) {
                double sum = 0.0;
                for (uint_t p : range(i)) {
                    sum += l.safe(i,p)*y.safe[p];
                }

                y.safe[i] = (x.safe[ipiv.safe[i]] - sum)/l.safe(i,i);
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
                        sum += u.safe(i,p)*z.safe[p];
                    }

                    z.safe[i] = (y.safe[i] - sum)/u.safe(i,i);
                }
            }

            // Apply pivot
            for (uint_t i : range(n)) {
                r.safe[ipiv.safe[i]] = z.safe[i];
            }

            return r;
        }
    };
}
}

#endif
