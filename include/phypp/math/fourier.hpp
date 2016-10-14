#ifndef PHYPP_MATH_FOURIER_HPP
#define PHYPP_MATH_FOURIER_HPP

#ifndef NO_FFTW
#include <fftw3.h>
#endif
#include "phypp/core/vec.hpp"
#include "phypp/math/complex.hpp"

namespace phypp {
    #ifndef NO_FFTW
    // Compute the Fast Fourier Transform (FFT) of the provided 2d array
    inline vec2cd fft(const vec2d& v) {
        vec2cd r(v.dims);

        fftw_plan p;
        p = fftw_plan_dft_r2c_2d(v.dims[0], v.dims[1],
            const_cast<double*>(v.data.data()),
            reinterpret_cast<fftw_complex*>(r.data.data()), FFTW_ESTIMATE);
        fftw_execute(p);
        fftw_destroy_plan(p);

        return r;
    }

    // Compute the Fast Fourier Transform (FFT) of the provided 2d array
    inline vec2d ifft(const vec2cd& v) {
        vec2d r(v.dims);

        fftw_plan p;
        p = fftw_plan_dft_c2r_2d(v.dims[0], v.dims[1],
            const_cast<fftw_complex*>(reinterpret_cast<const fftw_complex*>(v.data.data())),
            r.data.data(), FFTW_ESTIMATE);
        fftw_execute(p);
        fftw_destroy_plan(p);

        return r;
    }
    #endif
}

#endif
