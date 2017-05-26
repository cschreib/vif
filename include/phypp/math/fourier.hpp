#ifndef PHYPP_MATH_FOURIER_HPP
#define PHYPP_MATH_FOURIER_HPP

#ifndef NO_FFTW
#include <fftw3.h>
#endif
#include "phypp/core/vec.hpp"
#include "phypp/utility/thread.hpp"
#include "phypp/math/complex.hpp"

namespace phypp {
    namespace impl {
        inline std::mutex& fftw_planner_mutex() {
            static std::mutex m;
            return m;
        }
    }

    #ifndef NO_FFTW
    // Compute the Fast Fourier Transform (FFT) of the provided 2d array
    inline void fft(const vec2d& v, vec2cd& r) {
        fftw_plan p;
        {
            std::lock_guard<std::mutex> lock(impl::fftw_planner_mutex());
            p = fftw_plan_dft_r2c_2d(v.dims[0], v.dims[1],
                const_cast<double*>(v.data.data()),
                reinterpret_cast<fftw_complex*>(r.data.data()), FFTW_ESTIMATE);
        }

        fftw_execute(p);

        {
            std::lock_guard<std::mutex> lock(impl::fftw_planner_mutex());
            fftw_destroy_plan(p);
        }
    }

    // Compute the Fast Fourier Transform (FFT) of the provided 2d array
    inline vec2cd fft(const vec2d& v) {
        vec2cd r(v.dims);
        fft(v, r);
        return r;
    }

    // Compute the Fast Fourier Transform (FFT) of the provided 2d array
    // NB: the FFTW routine does not preserve the data in input, so the
    // input array has to be copied
    inline void ifft(vec2cd v, vec2d& r) {
        fftw_plan p;
        {
            std::lock_guard<std::mutex> lock(impl::fftw_planner_mutex());
            p = fftw_plan_dft_c2r_2d(v.dims[0], v.dims[1],
                reinterpret_cast<fftw_complex*>(v.data.data()),
                r.data.data(), FFTW_ESTIMATE);
        }

        fftw_execute(p);

        {
            std::lock_guard<std::mutex> lock(impl::fftw_planner_mutex());
            fftw_destroy_plan(p);
        }
    }

    // Compute the Fast Fourier Transform (FFT) of the provided 2d array
    // NB: the FFTW routine does not preserve the data in input, so the
    // input array has to be copied
    inline vec2d ifft(vec2cd v) {
        vec2d r(v.dims);
        ifft(v, r);
        return r;
    }
    #endif
}

#endif
