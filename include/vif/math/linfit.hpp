#ifndef VIF_MATH_LINFIT_HPP
#define VIF_MATH_LINFIT_HPP

#include "vif/core/vec.hpp"
#include "vif/core/error.hpp"
#include "vif/core/range.hpp"
#include "vif/math/base.hpp"
#include "vif/math/matrix.hpp"

namespace vif {
    struct linfit_result {
        bool          success = false;
        double        chi2 = dnan;
        vec1d         params;
        vec1d         errors;
        matrix::mat2d cov;

        // Reflection data
        MEMBERS1(success, chi2, params, errors, cov)
        MEMBERS2("linfit_result", MAKE_MEMBER(success), MAKE_MEMBER(chi2),
            MAKE_MEMBER(params), MAKE_MEMBER(errors), MAKE_MEMBER(cov))
    };

    namespace impl {
        template<typename T, typename TypeE>
        void linfit_make_cache_(vec2d& cache, const TypeE& ye, uint_t i, T&& t) {
            cache.safe(i,_) = flatten(t/ye);
        }

        template<typename T, typename TypeE, typename ... Args>
        void linfit_make_cache_(vec2d& cache, const TypeE& ye, uint_t i, T&& t, Args&& ... args) {
            cache.safe(i,_) = flatten(t/ye);
            linfit_make_cache_(cache, ye, i+1, args...);
        }

        template<typename TypeY, typename TypeE>
        linfit_result linfit_do_(const TypeY& y, const TypeE& ye, const vec2d& cache) {
            linfit_result fr;

            uint_t np = cache.dims[0];
            uint_t nm = cache.dims[1];

            // Solving 'y +/- e = sum over i of a[i]*x[i]' to get all a[i]'s
            matrix::mat2d alpha(np,np);
            vec1d beta(np);
            auto tmp = flatten(y/ye);
            for (uint_t i : range(np)) {
                for (uint_t j : range(np)) {
                    if (i <= j) {
                        alpha.safe(i,j) = 0.0;
                        // alpha(i,j) = sum over all points of x[i]*x[j]/e^2
                        for (uint_t m : range(nm)) {
                            alpha.safe(i,j) += cache.safe(i,m)*cache.safe(j,m);
                        }
                    } else {
                        alpha.safe(i,j) = alpha.safe(j,i);
                    }
                }

                beta.safe[i] = 0.0;
                // beta[i] = sum over all points of x[i]*y/e^2
                for (uint_t m : range(nm)) {
                    beta.safe[i] += cache.safe(i,m)*tmp.safe[m];
                }
            }

            if (!inplace_invert_symmetric(alpha)) {
                fr.success = false;
                fr.chi2 = dnan;
                fr.params = replicate(dnan, np);
                fr.errors = replicate(dnan, np);
                inplace_symmetrize(alpha);
                fr.cov = alpha;
                return fr;
            }

            inplace_symmetrize(alpha);
            fr.success = true;
            fr.params = alpha*beta;
            fr.errors = sqrt(diagonal(alpha));
            fr.cov = alpha;

            vec1d model(nm);
            for (uint_t m : range(nm)) {
                model.safe[m] = 0.0;
                for (uint_t i : range(np)) {
                    model.safe[m] += fr.params.safe[i]*cache.safe(i,m);
                }
            }

            fr.chi2 = total(sqr(model - tmp));

            return fr;
        }

        template<typename TY>
        void linfit_error_dims_(const TY& y, uint_t i) {}

        template<typename TY, typename U, typename ... Args>
        void linfit_error_dims_(const TY& y, uint_t i, const U& t, const Args& ... args) {
            vif_check(meta::same_dims_or_scalar(y, t), "incompatible dimensions between Y and X",
                i, " (", meta::dims(y), " vs. ", meta::dims(t), ")");
            linfit_error_dims_(y, i+1, args...);
        }
    }

    template<typename TypeY, typename TypeE, typename ... Args>
    linfit_result linfit(const TypeY& y, const TypeE& ye, Args&&... args) {
        bool bad = !meta::same_dims_or_scalar(y, ye, args...);
        if (bad) {
            vif_check(meta::same_dims_or_scalar(y, ye), "incompatible dimensions between Y and "
                "YE arrays (", meta::dims(y), " vs. ", meta::dims(ye), ")");
            impl::linfit_error_dims_(y, 0, args...);
        }

        uint_t np = sizeof...(Args);
        uint_t nm = meta::size(y);

        vec2d cache(np,nm);
        impl::linfit_make_cache_(cache, ye, 0, std::forward<Args>(args)...);

        return impl::linfit_do_(y, ye, cache);
    }

    template<std::size_t Dim, typename TypeY, typename TypeE, typename TypeX>
    linfit_result linfit_pack(const vec<Dim,TypeY>& y, const vec<Dim,TypeE>& ye,
        const vec<Dim+1,TypeX>& x) {
        bool good = true;
        for (uint_t i : range(Dim)) {
            if (x.dims[i+1] != ye.dims[i] || x.dims[i+1] != y.dims[i]) {
                good = false;
                break;
            }
        }

        vif_check(good, "incompatible dimensions between X, Y and YE arrays (", x.dims,
            " vs. ", y.dims, " vs. ", ye.dims, ")")

        linfit_result fr;

        uint_t np = x.dims[0];
        uint_t nm = y.size();

        vec2d cache(np,nm);
        for (uint_t i : range(np))
        for (uint_t j : range(nm)) {
            cache.safe(i,j) = x.safe[i*x.pitch(0) + j]/ye.safe[j];
        }

        return impl::linfit_do_(y, ye, cache);
    }

    template<typename TypeE>
    struct linfit_batch_t {
        TypeE         ye;
        vec2d         cache;
        vec1d         beta;
        matrix::mat2d alpha;
        linfit_result fr;

        struct pack_tag {};

    private :
        void update_matrix_() {
            uint_t np = cache.dims[0];
            uint_t nm = cache.dims[1];

            beta.resize(np);
            alpha.resize(np,np);

            // Solving 'y +/- e = sum over i of a[i]*x[i]' to get all a[i]'s
            for (uint_t i : range(np))
            for (uint_t j : range(np)) {
                if (i <= j) {
                    alpha.safe(i,j) = 0.0;
                    // alpha(i,j) = sum over all points of x[i]*x[j]/e^2
                    for (uint_t m : range(nm)) {
                        alpha.safe(i,j) += cache.safe(i,m)*cache.safe(j,m);
                    }
                } else {
                    alpha.safe(i,j) = alpha.safe(j,i);
                }
            }
        }

        void update_matrix_(uint_t i) {
            uint_t np = cache.dims[0];
            uint_t nm = cache.dims[1];

            for (uint_t j : range(np)) {
                alpha.safe(i,j) = 0.0;
                for (uint_t m : range(nm)) {
                    alpha.safe(i,j) += cache.safe(i,m)*cache.safe(j,m);
                }
                alpha.safe(j,i) = alpha.safe(i,j);
            }
        }

        void invert_matrix_() {
            uint_t np = cache.dims[0];
            uint_t nm = cache.dims[1];

            if (!invert_symmetric(alpha, fr.cov)) {
                fr.success = false;
                fr.chi2 = dnan;
                fr.params = replicate(dnan, np);
                fr.errors = replicate(dnan, np);
            } else {
                fr.success = true;
                fr.errors = sqrt(matrix::diagonal(fr.cov));
            }

            inplace_symmetrize(fr.cov);
        }

        template<typename TypeX>
        void update_model_(uint_t i, const TypeX& x) {
            uint_t np = cache.dims[0];
            uint_t nm = cache.dims[1];

            vif_check(i < np, "index goes beyond available models (", i, " vs. ", np, ")");

            for (uint_t j : range(nm)) {
                cache.safe(i,j) = x.safe[j]/ye.safe[j];
            }

            update_matrix_(i);
        }

    public:
        template<typename ... Args>
        linfit_batch_t(const TypeE& e, Args&&... args) : ye(e) {
            uint_t np = sizeof...(Args);
            uint_t nm = meta::size(ye);

            cache.resize(np,nm);
            impl::linfit_make_cache_(cache, ye, 0, std::forward<Args>(args)...);

            update_matrix_();
            invert_matrix_();
        }

        template<typename TypeX>
        linfit_batch_t(pack_tag, const TypeE& e, const TypeX& x) : ye(e) {
            uint_t np = x.dims[0];
            uint_t nm = meta::size(ye);

            cache.resize(np,nm);
            for (uint_t i : range(np))
            for (uint_t j : range(nm)) {
                cache.safe(i,j) = x.safe[i*x.pitch(0) + j]/ye.safe[j];
            }

            update_matrix_();
            invert_matrix_();
        }

        template<typename TypeX>
        void update_model(uint_t i, const TypeX& x) {
            update_model_(i, x);
            invert_matrix_();
        }

        void update_models() {
            invert_matrix_();
        }

        template<typename TypeX, typename ... Args>
        void update_models(uint_t i, const TypeX& x, const Args& ... args) {
            update_model_(i, x);
            update_models(args...);
        }

        template<typename TypeY>
        void fit_nochi2(const TypeY& y) {
            vif_check(meta::same_dims_or_scalar(y, ye), "incompatible dimensions between Y and "
                "YE arrays (", meta::dims(y), " vs. ", meta::dims(ye), ")");

            if (!fr.success) return;

            uint_t np = cache.dims[0];
            uint_t nm = cache.dims[1];

            auto tmp = flatten(y/ye);
            for (uint_t i : range(np)) {
                beta.safe[i] = 0.0;
                // beta[i] = sum over all points of x[i]*y/e^2
                for (uint_t m : range(nm)) {
                    beta.safe[i] += cache.safe(i,m)*tmp.safe[m];
                }
            }

            fr.params = fr.cov*beta;
        }

        template<typename TypeY>
        void fit(const TypeY& y) {
            fit_nochi2(y);
            update_chi2(y);
        }

        template<typename TypeY>
        void update_chi2(const TypeY& y) {
            uint_t np = cache.dims[0];
            uint_t nm = cache.dims[1];

            vec1d model(nm);
            for (uint_t m : range(nm)) {
                model.safe[m] = 0.0;
                for (uint_t i : range(np)) {
                    model.safe[m] += fr.params.safe[i]*cache.safe(i,m);
                }
            }

            fr.chi2 = total(sqr(model - flatten(y/ye)));
        }
    };

    template<typename TypeE, typename ... Args>
    linfit_batch_t<TypeE> linfit_batch(const TypeE& e, Args&&... args) {
        return linfit_batch_t<TypeE>(e, std::forward<Args>(args)...);
    }

    template<typename TypeE, typename ... Args>
    linfit_batch_t<TypeE> linfit_pack_batch(const TypeE& e, Args&&... args) {
        return linfit_batch_t<TypeE>(typename linfit_batch_t<TypeE>::pack_tag{},
            e, std::forward<Args>(args)...);
    }
}

#endif
