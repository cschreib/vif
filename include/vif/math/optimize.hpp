#ifndef VIF_MATH_OPTIMIZE_HPP
#define VIF_MATH_OPTIMIZE_HPP

#include "vif/core/vec.hpp"
#include "vif/core/error.hpp"
#include "vif/core/range.hpp"

#ifndef NO_GSL
#include <gsl/gsl_multimin.h>
#endif

namespace vif {
    struct minimize_params {
        double initial_step = 0.1;
        double ftol = 1e-3;
        uint_t max_iter = 1000;
    };

    struct minimize_result {
        bool success = false;
        vec1d params;
        double value = dnan;
        uint_t niter = 0;
    };

    enum class minimize_function_output {
        value, derivatives, all
    };

    template<typename F>
    minimize_result minimize_bfgs(const minimize_params& opts, const vec1d& start, F&& func) {
        using func_ptr = typename std::decay<decltype(func)>::type*;

        // Initialize return value
        minimize_result ret;
        ret.params.resize(start.dims);

        const uint_t n = start.size();

        // Initialize minimized function
        gsl_multimin_function_fdf mf;
        mf.n = n;
        mf.params = reinterpret_cast<void*>(&func);

        mf.f = [](const gsl_vector* x, void* p) {
            const uint_t tn = x->size;
            vec1d vv(tn);
            for (uint_t i : range(tn)) {
                vv.safe[i] = gsl_vector_get(x, i);
            }

            func_ptr fp = reinterpret_cast<func_ptr>(p);
            vec1d r = (*fp)(vv, minimize_function_output::value);

            return r.safe[0];
        };

        mf.df = [](const gsl_vector* x, void* p, gsl_vector* g) {
            const uint_t tn = x->size;
            vec1d vv(tn);
            for (uint_t i : range(tn)) {
                vv.safe[i] = gsl_vector_get(x, i);
            }

            func_ptr fp = reinterpret_cast<func_ptr>(p);
            vec1d r = (*fp)(vv, minimize_function_output::derivatives);

            for (uint_t i : range(tn)) {
                gsl_vector_set(g, i, r.safe[i+1]);
            }
        };

        mf.fdf = [](const gsl_vector* x, void* p, double* f, gsl_vector* g) {
            const uint_t tn = x->size;
            vec1d vv(tn);
            for (uint_t i : range(tn)) {
                vv.safe[i] = gsl_vector_get(x, i);
            }

            func_ptr fp = reinterpret_cast<func_ptr>(p);
            vec1d r = (*fp)(vv, minimize_function_output::all);

            *f = r.safe[0];

            for (uint_t i : range(tn)) {
                gsl_vector_set(g, i, r.safe[i+1]);
            }
        };

        gsl_vector* x = nullptr;
        gsl_multimin_fdfminimizer* m = nullptr;

        try {
            // Initial conditions
            x = gsl_vector_alloc(n);
            for (uint_t i : range(n)) {
                gsl_vector_set(x, i, start.safe[i]);
            }

            // Initialize minimizer
            const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_vector_bfgs2;
            m = gsl_multimin_fdfminimizer_alloc(T, n);
            gsl_multimin_fdfminimizer_set(m, &mf, x, opts.initial_step, opts.ftol);

            // Minimization loop
            int status;

            do {
                ++ret.niter;

                status = gsl_multimin_fdfminimizer_iterate(m);
                if (status != 0) break;

                status = gsl_multimin_test_gradient(m->gradient, 1e-3);
                if (status == GSL_SUCCESS) {
                    ret.success = true;
                    break;
                }
            } while (status == GSL_CONTINUE && ret.niter < opts.max_iter);

            // Write output
            for (uint_t i : range(n)) {
                ret.params.safe[i] = gsl_vector_get(m->x, i);
            }

            ret.value = m->f;

            // Cleanup
            gsl_multimin_fdfminimizer_free(m);
            gsl_vector_free(x);
        } catch (...) {
            // Cleanup
            if (m) gsl_multimin_fdfminimizer_free(m);
            if (x) gsl_vector_free(x);
        }

        return ret;
    }
}

#endif
