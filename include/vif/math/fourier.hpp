#ifndef VIF_MATH_FOURIER_HPP
#define VIF_MATH_FOURIER_HPP

#ifndef NO_FFTW
#include <fftw3.h>
#endif
#include "vif/core/vec.hpp"
#include "vif/utility/thread.hpp"
#include "vif/math/complex.hpp"

namespace vif {
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
                const_cast<double*>(v.raw_data()),
                reinterpret_cast<fftw_complex*>(r.raw_data()), FFTW_ESTIMATE);
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

    // Compute the Fast Fourier Transform (FFT) of the provided 2d complex array
    inline void fft_c2c(const vec2cd& v, vec2cd& r) {
        fftw_plan p;
        {
            std::lock_guard<std::mutex> lock(impl::fftw_planner_mutex());
            p = fftw_plan_dft_2d(v.dims[0], v.dims[1],
                const_cast<fftw_complex*>(reinterpret_cast<const fftw_complex*>(v.raw_data())),
                reinterpret_cast<fftw_complex*>(r.raw_data()), FFTW_FORWARD, FFTW_ESTIMATE);
        }

        fftw_execute(p);

        {
            std::lock_guard<std::mutex> lock(impl::fftw_planner_mutex());
            fftw_destroy_plan(p);
        }
    }

    // Compute the Fast Fourier Transform (FFT) of the provided 2d array
    inline vec2cd fft_c2c(const vec2cd& v) {
        vec2cd r(v.dims);
        fft_c2c(v, r);
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
                reinterpret_cast<fftw_complex*>(v.raw_data()),
                r.raw_data(), FFTW_ESTIMATE);
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

    // Compute the Fast Fourier Transform (FFT) of the provided 2d array
    // NB: the FFTW routine does not preserve the data in input, so the
    // input array has to be copied
    inline void ifft_c2c(vec2cd v, vec2cd& r) {
        fftw_plan p;
        {
            std::lock_guard<std::mutex> lock(impl::fftw_planner_mutex());
            p = fftw_plan_dft_2d(v.dims[0], v.dims[1],
                reinterpret_cast<fftw_complex*>(v.raw_data()),
                reinterpret_cast<fftw_complex*>(r.raw_data()), FFTW_BACKWARD, FFTW_ESTIMATE);
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
    inline vec2cd ifft_c2c(vec2cd v) {
        vec2cd r(v.dims);
        ifft_c2c(v, r);
        return r;
    }
    #endif
}

#endif
