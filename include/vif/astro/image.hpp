#ifndef VIF_ASTRO_IMAGE_HPP
#define VIF_ASTRO_IMAGE_HPP

#include "vif/core/vec.hpp"
#include "vif/core/error.hpp"
#include "vif/core/range.hpp"
#include "vif/utility/generic.hpp"
#include "vif/math/base.hpp"
#include "vif/math/fourier.hpp"
#include "vif/astro/wcs.hpp"

namespace vif {
namespace astro {
    template<typename Type>
    vec<2,meta::rtype_t<Type>> enlarge(const vec<2,Type>& v, const std::array<uint_t,4> upix,
        const meta::rtype_t<Type>& def = 0.0) {

        vec<2,meta::rtype_t<Type>> r = replicate(def,
            v.dims[0]+upix[0]+upix[2], v.dims[1]+upix[1]+upix[3]);

        for (uint_t x : range(v.dims[0]))
        for (uint_t y : range(v.dims[1])) {
            r.safe(x+upix[0],y+upix[1]) = v.safe(x,y);
        }

        return r;
    }

    template<typename Type>
    vec<2,meta::rtype_t<Type>> enlarge(const vec<2,Type>& v, uint_t upix,
        const meta::rtype_t<Type>& def = 0.0) {
        return enlarge(v, {{upix, upix, upix, upix}}, def);
    }

    template<typename Type>
    vec<2,meta::rtype_t<Type>> shrink(const vec<2,Type>& v, const std::array<uint_t,4> upix) {
        if (upix[0]+upix[2] >= v.dims[0] || upix[1]+upix[3] >= v.dims[1]) {
            return vec<2,meta::rtype_t<Type>>();
        }

        vec<2,meta::rtype_t<Type>> r(v.dims[0]-upix[0]-upix[2], v.dims[1]-upix[1]-upix[3]);

        for (uint_t x : range(r.dims[0]))
        for (uint_t y : range(r.dims[1])) {
            r.safe(x,y) = v.safe(x+upix[0],y+upix[1]);
        }

        return r;
    }

    template<typename Type>
    vec<2,meta::rtype_t<Type>> shrink(const vec<2,Type>& v, uint_t upix) {
        return shrink(v, {{upix, upix, upix, upix}});
    }

    // Get a sub region 'reg' inside an image 'v', with reg = {x0, y0, x1, y1} (inclusive)
    // assuming that "x" and "y" correspond to v(x,y).
    // The sub region is returned as two index vectors, the first containing indices in the
    // image 'v', and the second containing indices in the region 'reg'. The number of
    // returned indices may be smaller than the size of the requested region, if a fraction of
    // the region falls out of the image boundaries. In particular, the index vectors will be
    // empty if there is no overlap between the image and the requested region.
    template<typename TypeV, typename TypeR = int_t>
    void subregion(const vec<2,TypeV>& v, const vec<1,TypeR>& reg, vec1u& rv, vec1u& rr) {
        vif_check(reg.size() == 4, "invalid region parameter "
            "(expected 4 components, got ", reg.size(), ")");

        int_t nvx = v.dims[0], nvy = v.dims[1];
        int_t nx = reg.safe[2]-reg.safe[0]+1, ny = reg.safe[3]-reg.safe[1]+1;

        vec1i vreg = reg;
        vec1i sreg = {0,0,nx-1,ny-1};

        // Early exit when not covered
        if (reg.safe[0] >= nvx || reg.safe[2] < 0 || reg.safe[0] > reg.safe[2] ||
            reg.safe[1] >= nvy || reg.safe[3] < 0 || reg.safe[1] > reg.safe[3]) {
            rv.clear(); rr.clear();
            return;
        }

        if (reg.safe[0] < 0) {
            vreg.safe[0] = 0;
            sreg.safe[0] += 0 - reg.safe[0];
        }
        if (reg.safe[2] >= nvx) {
            vreg.safe[2] = nvx-1;
            sreg.safe[2] -= reg.safe[2] - (nvx-1);
        }

        if (reg.safe[1] < 0) {
            vreg.safe[1] = 0;
            sreg.safe[1] += 0 - reg.safe[1];
        }
        if (reg.safe[3] >= nvy) {
            vreg.safe[3] = nvy-1;
            sreg.safe[3] -= reg.safe[3] - (nvy-1);
        }

        int_t nnx = vreg.safe[2]-vreg.safe[0]+1;
        int_t nny = vreg.safe[3]-vreg.safe[1]+1;
        uint_t npix = nnx*nny;

        rv.resize(npix);
        rr.resize(npix);
        for (uint_t i : range(npix)) {
            rv.safe[i] = (i%nny) + vreg.safe[1] + (i/nny + vreg.safe[0])*nvy;
            rr.safe[i] = (i%nny) + sreg.safe[1] + (i/nny + sreg.safe[0])*ny;
        }
    }

    template<typename TypeV, typename TypeR = int_t>
    typename vec<2,TypeV>::effective_type subregion(const vec<2,TypeV>& v,
        const vec<1,TypeR>& reg, const typename vec<2,TypeV>::rtype& def = 0.0) {

        vec1u rr, rs;
        subregion(v, reg, rr, rs);

        int_t nx = reg(2)-reg(0)+1, ny = reg(3)-reg(1)+1;
        vec<2,meta::rtype_t<TypeV>> sub = replicate(meta::rtype_t<TypeV>(def), nx, ny);

        sub.safe[rs] = v.safe[rr];

        return sub;
    }

    template<typename TypeV>
    typename vec<2,TypeV>::effective_type translate(const vec<2,TypeV>& v, double dx, double dy,
        typename vec<2,TypeV>::rtype def = 0.0) {

        vec<2,meta::rtype_t<TypeV>> trs(v.dims);
        for (uint_t x : range(v.dims[0]))
        for (uint_t y : range(v.dims[1])) {
            trs.safe(x,y) = bilinear_strict(v, x - dx, y - dy, def);
        }

        return trs;
    }

    template<typename TypeV>
    typename vec<2,TypeV>::effective_type translate_bicubic(const vec<2,TypeV>& v, double dx, double dy,
        typename vec<2,TypeV>::rtype def = 0.0) {

        vec<2,meta::rtype_t<TypeV>> trs = v;
        vec1d t = indgen<double>(v.dims[1]) - dy;
        for (uint_t x : range(v.dims[0])) {
            trs.safe(x,_) = interpolate_3spline(trs.safe(x,_), t);
        }

        t = indgen<double>(v.dims[0]) - dx;
        for (uint_t y : range(v.dims[1])) {
            trs.safe(_,y) = interpolate_3spline(trs.safe(_,y), t);
        }

        return trs;
    }

    template<typename TypeV>
    typename vec<2,TypeV>::effective_type flip_x(const vec<2,TypeV>& v) {
        auto r = v.concretise();
        for (uint_t y : range(v.dims[0]))
        for (uint_t x : range(v.dims[1]/2)) {
            std::swap(r.safe(y,x), r.safe(y,v.dims[1]-1-x));
        }

        return r;
    }

    template<typename TypeV>
    typename vec<2,TypeV>::effective_type flip_y(const vec<2,TypeV>& v) {
        auto r = v.concretise();
        for (uint_t y : range(v.dims[0]/2))
        for (uint_t x : range(v.dims[1])) {
            std::swap(r.safe(y,x), r.safe(v.dims[0]-1-y,x));
        }

        return r;
    }

    template<typename TypeV, typename TypeD = double>
    typename vec<2,TypeV>::effective_type scale(const vec<2,TypeV>& v, double factor,
        typename vec<2,TypeV>::rtype def = 0.0) {

        auto r = v.concretise();

        for (int_t y : range(v.dims[0]))
        for (int_t x : range(v.dims[1])) {
            r.safe(uint_t(y),uint_t(x)) = bilinear_strict(v,
                (y-int_t(v.dims[0]/2))/factor + v.dims[0]/2,
                (x-int_t(v.dims[1]/2))/factor + v.dims[1]/2,
                def
            );
        }

        return r;
    }

    template<typename TypeV, typename TypeD = double>
    typename vec<2,TypeV>::effective_type scale_bicubic(const vec<2,TypeV>& v, double factor,
        typename vec<2,TypeV>::rtype def = 0.0) {

        auto r = v.concretise();

        for (int_t y : range(v.dims[0]))
        for (int_t x : range(v.dims[1])) {
            r.safe(uint_t(y),uint_t(x)) = bicubic_strict(v,
                (y-int_t(v.dims[0]/2))/factor + v.dims[0]/2,
                (x-int_t(v.dims[1]/2))/factor + v.dims[1]/2,
                def
            );
        }

        return r;
    }

    template<typename TypeV, typename TypeD = double>
    typename vec<2,TypeV>::effective_type rotate(const vec<2,TypeV>& v, double angle,
        TypeD def = 0.0) {

        auto r = v.concretise();

        double ca = cos(angle*dpi/180.0);
        double sa = sin(angle*dpi/180.0);

        for (int_t y : range(v.dims[0]))
        for (int_t x : range(v.dims[1])) {
            double dy = (y-int_t(v.dims[0]/2));
            double dx = (x-int_t(v.dims[1]/2));
            r.safe(uint_t(y),uint_t(x)) = bilinear_strict(v,
                dy*ca - dx*sa + v.dims[0]/2,
                dx*ca + dy*sa + v.dims[1]/2,
                def
            );
        }

        return r;
    }

    template<typename TypeV, typename TypeD = double>
    typename vec<2,TypeV>::effective_type rotate_bicubic(const vec<2,TypeV>& v, double angle,
        TypeD def = 0.0) {

        auto r = v.concretise();

        double ca = cos(angle*dpi/180.0);
        double sa = sin(angle*dpi/180.0);

        for (int_t y : range(v.dims[0]))
        for (int_t x : range(v.dims[1])) {
            double dy = (y-int_t(v.dims[0]/2));
            double dx = (x-int_t(v.dims[1]/2));
            r.safe(uint_t(y),uint_t(x)) = bicubic_strict(v,
                dy*ca - dx*sa + v.dims[0]/2,
                dx*ca + dy*sa + v.dims[1]/2,
                def
            );
        }

        return r;
    }

    template<typename T>
    void inplace_shift(vec<2,T>& v, int_t tsx, int_t tsy) {
        tsx = (-tsx) % int_t(v.dims[0]);
        tsy = (-tsy) % int_t(v.dims[1]);
        if (tsx < 0) tsx = int_t(v.dims[0])+tsx;
        if (tsy < 0) tsy = int_t(v.dims[1])+tsy;

        if (tsy != 0) {
            for (uint_t ix : range(v.dims[0])) {
                auto st = v.raw_stride(ix,_);
                std::rotate(st.begin(), st.begin() + tsy, st.end());
            }
        }

        if (tsx != 0) {
            for (uint_t iy : range(v.dims[1])) {
                auto st = v.raw_stride(_,iy);
                std::rotate(st.begin(), st.begin() + tsx, st.end());
            }
        }
    }

    template<typename T>
    vec<2,T> shift(vec<2,T> v, int_t tsx, int_t tsy) {
        inplace_shift(v, tsx, tsy);
        return v;
    }

    template<typename T>
    void inplace_translate_integer(vec<2,T>& v, int_t tsx, int_t tsy, T def = 0) {
        static_assert(!std::is_pointer<T>::value, "this algorithm cannot work with views");

        tsx = (-tsx) % int_t(v.dims[0]);
        tsy = (-tsy) % int_t(v.dims[1]);

        if (tsy > 0) {
            for (uint_t ix : range(v.dims[0])) {
                auto st = v.raw_stride(ix,_);
                auto it = std::copy(st.begin() + tsy, st.end(), st.begin());
                std::fill(it, st.end(), def);
            }
        } else if (tsy < 0) {
            for (uint_t ix : range(v.dims[0])) {
                auto st = v.raw_stride(ix,_);
                auto it = std::copy_backward(st.begin(), st.end() + tsy, st.end());
                std::fill(st.begin(), it, def);
            }
        }

        if (tsx > 0) {
            for (uint_t iy : range(v.dims[1])) {
                auto st = v.raw_stride(_,iy);
                auto it = std::copy(st.begin() + tsx, st.end(), st.begin());
                std::fill(it, st.end(), def);
            }
        } else if (tsx < 0) {
            for (uint_t iy : range(v.dims[1])) {
                auto st = v.raw_stride(_,iy);
                auto it = std::copy_backward(st.begin(), st.end() + tsx, st.end());
                std::fill(st.begin(), it, def);
            }
        }
    }

    template<typename T>
    vec<2,T> translate_integer(vec<2,T> v, int_t tsx, int_t tsy, T def = 0) {
        inplace_translate_integer(v, tsx, tsy, def);
        return v;
    }

    template<typename T>
    vec<2,meta::rtype_t<T*>> translate_integer(vec<2,T*> v, int_t tsx, int_t tsy, T def = 0) {
        vec<2,meta::rtype_t<T*>> tv = v.concretise();
        inplace_translate_integer(tv, tsx, tsy, def);
        return tv;
    }

    template<typename T>
    vec<2,meta::rtype_t<T>> recenter(const vec<2,T>& img, int_t cy, int_t cx, std::array<uint_t,2> dims,
        meta::rtype_t<T> def = 0) {

        vec<2,meta::rtype_t<T>> nimg = replicate(def, dims);
        vec1u ip, it;
        int_t ty0 = cy - int_t(dims[0])/2;
        int_t tx0 = cx - int_t(dims[1])/2;
        int_t ty1 = ty0 + dims[0]-1;
        int_t tx1 = tx0 + dims[1]-1;
        subregion(img, {ty0, tx0, ty1, tx1}, ip, it);
        nimg[it] = img[ip];
        return nimg;
    }

    template<typename T>
    vec<2,meta::rtype_t<T>> recenter(const vec<2,T>& img, int_t cy, int_t cx, meta::rtype_t<T> def = 0) {
        return translate_integer(img, int_t(img.dims[0])/2 - cy, int_t(img.dims[1])/2 - cx, def);
    }

    inline vec2d circular_mask(const std::array<uint_t,2>& dims, double radius, double x, double y) {
        vec2d m(dims);

        vif_check(radius >= 0, "radius must be a positive number");

        // Identify the needed region
        if (x+radius > 0 && y+radius > 0) {
            uint_t x0 = floor(x - radius) > 0         ? floor(x - radius) : 0;
            uint_t y0 = floor(y - radius) > 0         ? floor(y - radius) : 0;
            uint_t x1 = ceil(x + radius)  < dims[0]-1 ? ceil(x + radius)  : dims[0]-1;
            uint_t y1 = ceil(y + radius)  < dims[1]-1 ? ceil(y + radius)  : dims[1]-1;

            radius *= radius;
            for (uint_t ix = x0; ix <= x1; ++ix)
            for (uint_t iy = y0; iy <= y1; ++iy) {
                m.safe(ix,iy) = 0.25*(
                    (sqr(ix-0.5 - x) + sqr(iy-0.5 - y) <= radius) +
                    (sqr(ix+0.5 - x) + sqr(iy-0.5 - y) <= radius) +
                    (sqr(ix+0.5 - x) + sqr(iy+0.5 - y) <= radius) +
                    (sqr(ix-0.5 - x) + sqr(iy+0.5 - y) <= radius)
                );
            }
        }

        return m;
    }

    inline vec2d circular_mask(const std::array<uint_t,2>& dims, double radius) {
        return circular_mask(dims, radius, dims[0]/2.0, dims[1]/2.0);
    }

    template<typename Type>
    vec<1, meta::rtype_t<Type>> radial_profile(const vec<2,Type>& img, uint_t npix) {
        // TODO: (optimization) really...
        vec<1, meta::rtype_t<Type>> res(npix);
        uint_t hsx = img.dims[0]/2;
        uint_t hsy = img.dims[1]/2;
        res[0] = img(hsx,hsy);
        for (uint_t i : range(1u, npix)) {
            vec2d mask = circular_mask(img.dims, double(hsx), double(hsy), i)*
                (1.0 - circular_mask(img.dims, double(hsx), double(hsy), i-1));
            res.safe[i] = total(mask*img)/total(mask);
        }

        return res;
    }

    template<typename F>
    auto generate_img(const std::array<uint_t,2>& dims, F&& expr) -> vec<2,decltype(expr(0,0))> {
        vec<2,decltype(expr(0,0))> img(dims);
        for (uint_t x : range(img.dims[0]))
        for (uint_t y : range(img.dims[1])) {
            img.safe(x,y) = expr(x,y);
        }

        return img;
    }

    inline vec2d gaussian_profile(const std::array<uint_t,2>& dims, double sigma,
        double x0, double y0) {

        return generate_img(dims, [&](double x, double y) {
            return integrate_gauss_2d(x-0.5, x+0.5, y-0.5, y+0.5, x0, y0, sigma);
        });
    }

    inline vec2d gaussian_profile(const std::array<uint_t,2>& dims, double sigma) {
        return gaussian_profile(dims, sigma, dims[0]/2, dims[1]/2);
    }

    // Perform the convolution of two 2D arrays, assuming the second one is the kernel.
    // Naive loop implementation: this is slow but simple and reliable.
    template<typename TypeY1, typename TypeY2>
    auto convolve2d_naive(const vec<2,TypeY1>& map, const vec<2,TypeY2>& kernel) ->
        vec<2,decltype(map[0]*kernel[0])> {
        vif_check(kernel.dims[0]%2 == 1 && kernel.dims[1]%2 == 1,
            "kernel must have odd dimensions (", kernel.dims, ")");

        uint_t hxsize = kernel.dims[0]/2;
        uint_t hysize = kernel.dims[1]/2;

        vec<2,decltype(map[0]*kernel[0])> r(map.dims);
        for (uint_t kx = 0; kx < kernel.dims[0]; ++kx)
        for (uint_t ky = 0; ky < kernel.dims[1]; ++ky) {
            const auto kw = kernel.safe(kx, ky);
            if (kw == 0.0) continue;

            uint_t x0 = (kx >= hxsize ? 0 : hxsize-kx);
            uint_t xn = map.dims[0] - (kx >= hxsize ? kx-hxsize : 0);
            uint_t y0 = (ky >= hysize ? 0 : hysize-ky);
            uint_t yn = map.dims[1] - (ky >= hysize ? ky-hysize : 0);

            for (uint_t x = x0; x < xn; ++x)
            for (uint_t y = y0; y < yn; ++y) {
                r.safe(x+kx-hxsize,y+ky-hysize) += map.safe(x,y)*kw;
            }
        }

        return r;
    }

    // Perform the convolution of two 2D arrays, assuming the second one is the kernel.
    // Note: If the FFTW library is not used, falls back to convolve2d_naive().
    template<typename TypeY1, typename TypeY2>
    vec2cd convolve2d_prepare_kernel(const vec<2,TypeY1>& map, const vec<2,TypeY2>& kernel,
        uint_t& hsx, uint_t& hsy) {
#ifdef NO_FFTW
        static_assert(!std::is_same<TypeY1,TypeY1>::value, "this function requires the FFTW "
            "library");
        return vec2cd();
#else
        vif_check(kernel.dims[0]%2 == 1 && kernel.dims[1]%2 == 1,
            "kernel must have odd dimensions (", kernel.dims, ")");

        hsx = kernel.dims[0]/2; hsy = kernel.dims[1]/2;

        // Resize kernel to map size, with kernel center at (0,0), and some padding
        // to avoid cyclic borders (assume image is 0 outside)
        vec2d tkernel = enlarge(kernel, {{0, 0, map.dims[0]-1, map.dims[1]-1}});
        inplace_shift(tkernel, -hsx, -hsy);

        // Put the kernel in Fourier space
        return fft(tkernel);
#endif
    }

    // Perform the convolution of a 2D array with a kernel (given in Fourier space)
    template<typename TypeY1>
    auto convolve2d(const vec<2,TypeY1>& map, const vec2cd& kernel, uint_t hsx, uint_t hsy) ->
        vec<2,meta::rtype_t<TypeY1>> {
#ifdef NO_FFTW
        static_assert(!std::is_same<TypeY1,TypeY1>::value, "this function requires the FFTW "
            "library");
        return map;
#else
        // Pad image to prevent issues with cyclic borders
        vec2d tmap = enlarge(map, {{hsx, hsy, hsx, hsy}});

        // Perform the convolution in Fourier space
        auto cimg = fft(tmap)*kernel;

        // Go back to real space and shrink map back to original dimensions
        return shrink(ifft(cimg), {{hsx, hsy, hsx, hsy}})/cimg.size();
#endif
    }

    // Perform the convolution of two 2D arrays, assuming the second one is the kernel.
    // Note: If the FFTW library is not used, falls back to convolve2d_naive().
    template<typename TypeY1, typename TypeY2>
    auto convolve2d(const vec<2,TypeY1>& map, const vec<2,TypeY2>& kernel) ->
        vec<2,decltype(map[0]*kernel[0])> {
#ifdef NO_FFTW
        return convolve2d_naive(map, kernel);
#else
        vif_check(kernel.dims[0]%2 == 1 && kernel.dims[1]%2 == 1,
            "kernel must have odd dimensions (", kernel.dims, ")");

        uint_t hsx = kernel.dims[0]/2, hsy = kernel.dims[1]/2;

        // Pad image to prevent issues with cyclic borders
        vec2d tmap = enlarge(map, {{hsx, hsy, hsx, hsy}});

        // Resize kernel to map size, with kernel center at (0,0)
        vec2d tkernel(tmap.dims);

        // TODO: optimize this and catch case where kernel is larger than image
        vec1u px1 = hsx + indgen<uint_t>(hsx);
        vec1u py1 = hsy + indgen<uint_t>(hsy);
        vec1u px2 = hsx - 1 - indgen<uint_t>(hsx);
        vec1u py2 = hsy - 1 - indgen<uint_t>(hsy);

        vec1u ix1 = indgen<uint_t>(hsx);
        vec1u iy1 = indgen<uint_t>(hsy);
        vec1u ix2 = tmap.dims[0] - 1 - indgen<uint_t>(hsx);
        vec1u iy2 = tmap.dims[1] - 1 - indgen<uint_t>(hsy);

        tkernel(ix1,iy1) = kernel(px1,py1);
        tkernel(ix2,iy1) = kernel(px2,py1);
        tkernel(ix1,iy2) = kernel(px1,py2);
        tkernel(ix2,iy2) = kernel(px2,py2);

        // Perform the convolution in Fourier space
        auto cimg = fft(tmap)*fft(tkernel);

        // Go back to real space and shrink map back to original dimensions
        return shrink(ifft(cimg), {{hsx, hsy, hsx, hsy}})/cimg.size();
#endif
    }

#ifndef NO_FFTW
    struct convolver2d {
        const vec2d& kernel_normal;
        vec2cd kernel_fourier;
        uint_t hsx = 0, hsy = 0;
        fftw_plan pf, pi;
        bool pf_built = false, pi_built = false;
        vec2d tmap;
        vec2cd cimg;

        explicit convolver2d(const vec2d& k) : kernel_normal(k) {
            vif_check(k.dims[0]%2 == 1 && k.dims[1]%2 == 1,
                "kernel must have odd dimensions (", k.dims, ")");
        }

        convolver2d(const convolver2d&) = delete;
        convolver2d& operator=(const convolver2d&) = delete;
        convolver2d(convolver2d&&) = default;
        convolver2d& operator=(convolver2d&&) = default;

        ~convolver2d() {
            std::lock_guard<std::mutex> lock(impl::fftw_planner_mutex());
            if (pf_built) {
                fftw_destroy_plan(pf);
            }
            if (pi_built) {
                fftw_destroy_plan(pi);
            }
        }

    private :
        void make_plan(const vec2d& v, vec2cd& r) {
            if (!pf_built) {
                std::lock_guard<std::mutex> lock(impl::fftw_planner_mutex());
                pf = fftw_plan_dft_r2c_2d(v.dims[0], v.dims[1],
                    const_cast<double*>(v.raw_data()),
                    reinterpret_cast<fftw_complex*>(r.raw_data()), FFTW_ESTIMATE);
                pf_built = true;
            }
        }

        void make_plan(vec2cd& v, vec2d& r) {
            if (!pi_built) {
                std::lock_guard<std::mutex> lock(impl::fftw_planner_mutex());
                pi = fftw_plan_dft_c2r_2d(v.dims[0], v.dims[1],
                    reinterpret_cast<fftw_complex*>(v.raw_data()),
                    r.raw_data(), FFTW_ESTIMATE);
                pi_built = true;
            }
        }

        void fft(const vec2d& v, vec2cd& r) {
            make_plan(v, r);

            fftw_execute_dft_r2c(pf,
                const_cast<double*>(v.raw_data()),
                reinterpret_cast<fftw_complex*>(r.raw_data())
            );
        }

        vec2cd fft(const vec2d& v) {
            vec2cd r(v.dims);
            fft(v, r);
            return r;
        }

        // Compute the Fast Fourier Transform (FFT) of the provided 2d array
        void ifft(vec2cd& v, vec2d& r) {
            make_plan(v, r);

            fftw_execute_dft_c2r(pi,
                reinterpret_cast<fftw_complex*>(v.raw_data()),
                r.raw_data()
            );
        }

        // Compute the Fast Fourier Transform (FFT) of the provided 2d array
        vec2d ifft(vec2cd& v) {
            vec2d r(v.dims);
            ifft(v, r);
            return r;
        }

    public :
        void inplace_convolve(vec2d& map) {
            if (kernel_fourier.empty()) {
                hsx = kernel_normal.dims[0]/2; hsy = kernel_normal.dims[1]/2;

                // Resize kernel to map size, with kernel center at (0,0), and some padding
                // to avoid cyclic borders (assume image is 0 outside)
                vec2d tkernel = enlarge(kernel_normal, {{0, 0, map.dims[0]-1, map.dims[1]-1}});
                inplace_shift(tkernel, -hsx, -hsy);

                // Put the kernel in Fourier space
                kernel_fourier = this->fft(tkernel);

                tmap.resize(kernel_fourier.dims);
                cimg.resize(kernel_fourier.dims);
            }

            // Pad image to prevent issues with cyclic borders
            for (uint_t ix : range(map.dims[0]))
            for (uint_t iy : range(map.dims[1])) {
                tmap.safe(hsx+ix,hsy+iy) = map.safe(ix,iy);
            }
            for (uint_t ix : range(hsx)) {
                for (uint_t iy : range(tmap.dims[0])) {
                    tmap.safe(ix,iy) = 0;
                }
                for (uint_t iy : range(tmap.dims[0])) {
                    tmap.safe(tmap.dims[0]-hsx+ix,iy) = 0;
                }
            }
            for (uint_t ix : range(map.dims[0])) {
                for (uint_t iy : range(hsy)) {
                    tmap.safe(hsx+ix,iy) = 0;
                }
                for (uint_t iy : range(hsy)) {
                    tmap.safe(hsx+ix,tmap.dims[1]-hsy+iy) = 0;
                }
            }

            // Perform the convolution in Fourier space
            this->fft(tmap, cimg);
            cimg *= kernel_fourier;

            // Go back to real space and shrink map back to original dimensions
            this->ifft(cimg, tmap);

            for (uint_t ix : range(map.dims[0]))
            for (uint_t iy : range(map.dims[1])) {
                map.safe(ix,iy) = tmap.safe(hsx+ix,hsy+iy)/cimg.size();
            }
        }

        vec2d convolve(vec2d map) {
            inplace_convolve(map);
            return map;
        }
    };
#else
    struct convolver2d {
        explicit convolver2d(const vec2d&) {}

        template<typename Dummy>
        vec2d convolve(const vec2d& r) {
            static_assert(!std::is_same<Dummy,Dummy>::value, "this function requires the FFTW "
                "library");

            return r;
        }
    };
#endif

    inline convolver2d batch_convolve2d(const vec2d& k) {
        return convolver2d(k);
    }

    // Perform the convolution of two 2D arrays, assuming the second one is the kernel.
    // Note: If the FFTW library is not used, falls back to convolve2d_naive().
    template<typename T = void>
    vec2d make_transition_kernel(const vec2d& from, const vec2d& to) {
#ifdef NO_FFTW
        static_assert(!std::is_same<T,T>::value, "this function needs the FFTW library to work");
        return vec2d();
#else
        vec2cd cfrom = fft(from);
        vec2cd cto = fft(to);
        cto(_-(cfrom.dims[1]/2),_) /= cfrom(_-(cfrom.dims[1]/2),_);
        return shift(ifft(cto), cfrom.dims[0]/2, cfrom.dims[1]/2);
#endif
    }

    template<typename T, typename F>
    auto boxcar(const vec<2,T>& img, uint_t hsize, F&& func) ->
        vec<2,decltype(func(flatten(img)))> {

        vec<2,decltype(func(flatten(img)))> res(img.dims);
        for (uint_t x = 0; x < img.dims[0]; ++x)
        for (uint_t y = 0; y < img.dims[1]; ++y) {
            uint_t x0 = x >= hsize ? x - hsize : 0;
            uint_t x1 = x < img.dims[0]-hsize ? x + hsize : img.dims[0]-1;
            uint_t y0 = y >= hsize ? y - hsize : 0;
            uint_t y1 = y < img.dims[1]-hsize ? y + hsize : img.dims[1]-1;

            vec<1,meta::rtype_t<T>> tmp;
            tmp.reserve((x1-x0)*(y1-y0));
            for (uint_t tx = x0; tx <= x1; ++tx)
            for (uint_t ty = y0; ty <= y1; ++ty) {
                tmp.push_back(img.safe(tx,ty));
            }

            res.safe(x,y) = func(tmp);
        }

        return res;
    }

    inline vec2b mask_inflate(vec2b m, uint_t d) {
        if (d != 0) {
            // Locate edges
            vec1u ids; ids.reserve(3*sqrt(m.size()));
            for (uint_t x : range(m.dims[0]))
            for (uint_t y : range(m.dims[1])) {
                if (!m.safe(x,y)) continue;

                if ((x > 0 && !m.safe(x-1,y)) || (x < m.dims[0]-1 && !m.safe(x+1,y)) ||
                    (y > 0 && !m.safe(x,y-1)) || (y < m.dims[1]-1 && !m.safe(x,y+1))) {
                    ids.push_back(flat_id(m, x, y));
                }
            }

            if (!ids.empty()) {
                // Make kernel
                // (will contain exactly 4*d*(d+1)/2 elements)
                vec1i x; x.reserve(2*d*(d+1));
                vec1i y; y.reserve(2*d*(d+1));

                for (int_t k : range(1, d+1))
                for (int_t l : range(4*k)) {
                    int_t ml = l%k;
                    if (l < k) {
                        x.push_back(k-ml);
                        y.push_back(-ml);
                    } else if (l < 2*k) {
                        x.push_back(-ml);
                        y.push_back(-k+ml);
                    } else if (l < 3*k) {
                        x.push_back(-k+ml);
                        y.push_back(ml);
                    } else {
                        x.push_back(ml);
                        y.push_back(k-ml);
                    }
                }

                // Apply kernel to edge points
                for (uint_t i : range(ids)) {
                    vec1u j = mult_ids(m, ids[i]);

                    for (int_t k : range(x)) {
                        int_t tx = j.safe[0] + x.safe[k];
                        int_t ty = j.safe[1] + y.safe[k];

                        if (tx >= 0 && uint_t(tx) < m.dims[0] && ty >= 0 && uint_t(ty) < m.dims[1]) {
                            m.safe(uint_t(tx), uint_t(ty)) = true;
                        }
                    }
                }
            }
        }

        return m;
    }

    template <typename F>
    void foreach_segment(const vec2u& seg, const vec1u& mids, F&& func) {
        vec2u visited(seg.dims);

        vec1u sids;
        vec1u ox, oy;
        uint_t curseg = 0;

        auto fetch_neighbors = [&](uint_t ty, uint_t tx) {
            auto check_add = [&](uint_t tty, uint_t ttx) {
                if (visited(tty,ttx) != curseg) {
                    visited(tty,ttx) = curseg;
                    ox.push_back(ttx);
                    oy.push_back(tty);
                }
            };

            if (ty != 0)             check_add(ty-1,tx);
            if (ty != seg.dims[0]-1) check_add(ty+1,tx);
            if (tx != 0)             check_add(ty,tx-1);
            if (tx != seg.dims[1]-1) check_add(ty,tx+1);
        };

        for (uint_t s : range(mids)) {
            sids.clear();
            oy.clear(); ox.clear();

            curseg = seg[mids[s]];

            vec1u tmid = mult_ids(seg, mids[s]);
            oy.push_back(tmid[0]);
            ox.push_back(tmid[1]);

            while (!ox.empty()) {
                vec1u tox = std::move(ox), toy = std::move(oy);
                ox.clear(); oy.clear();
                for (uint_t i : range(tox)) {
                    if (seg(toy[i],tox[i]) == curseg) {
                        sids.push_back(flat_id(seg, toy[i], tox[i]));
                        fetch_neighbors(toy[i], tox[i]);
                    }
                }
            }

            func(s, sids);
        }
    }

    struct segment_params {
        // Minimum number of pixels per segment
        uint_t min_area = 0u;
        // First ID used to place segments on the map
        uint_t first_id = 1u;
    };

    struct segment_output {
        // ID of the segment as given on the segmentation map
        vec1u id;
        // Size of a segment in number of pixels
        vec1u area;
        // Center
        vec1d px, py;
        // Flat index of the first value of a segment
        vec1u origin;
    };

    // Function to segment a binary or integer map into multiple contiguous components.
    // Does no de-blending, use segment_deblend if you need it. Values of 0 in the
    // input binary map are also 0 in the segmentation map.
    template <typename T, typename enable = typename std::enable_if<!std::is_pointer<T>::value>::type>
    vec2u segment(vec<2,T> map, segment_output& out, const segment_params& params = segment_params()) {
        vec2u smap(map.dims);

        vif_check(params.first_id > 0, "first ID must be > 0");

        uint_t id = params.first_id;

        std::vector<uint_t> oy, ox;
        for (uint_t y : range(map.dims[0]))
        for (uint_t x : range(map.dims[1])) {
            if (map.safe(y,x) == 0) continue;

            // Found a guy, create entry in the output
            out.id.push_back(id);
            out.area.push_back(1u);
            out.px.push_back(x);
            out.py.push_back(y);
            out.origin.push_back(flat_id(map, y, x));

            // Use an A* - like algorithm to navigate around and
            // figure out its extents
            oy.clear(); ox.clear();

            // This function sets a point as belonging to the current segment
            // and adds the neighboring points to a search list for inspection
            auto process_point = [&ox,&oy,&map,&smap,&out,id](uint_t ty, uint_t tx) {
                map.safe(ty,tx) = 0;
                smap.safe(ty,tx) = id;

                auto check_add = [&ox,&oy,&map,&out](uint_t tty, uint_t ttx) {
                    if (map.safe(tty,ttx) != 0) {
                        out.area.back() += 1u;
                        out.px.back() += ttx;
                        out.py.back() += tty;

                        oy.push_back(tty);
                        ox.push_back(ttx);
                    }
                };

                if (ty != 0)             check_add(ty-1,tx);
                if (ty != map.dims[0]-1) check_add(ty+1,tx);
                if (tx != 0)             check_add(ty,tx-1);
                if (tx != map.dims[1]-1) check_add(ty,tx+1);
            };

            // Add the first point to the segment
            process_point(y, x);

            // While there are points to inspect, process them as well
            while (!ox.empty()) {
                uint_t ty = oy.back(); oy.pop_back();
                uint_t tx = ox.back(); ox.pop_back();
                process_point(ty, tx);
            }

            // Average position
            out.px.back() /= out.area.back();
            out.py.back() /= out.area.back();

            ++id;
        }

        // Erase too small regions
        if (count(out.area < params.min_area) != 0) {
            vec1u eids;
            foreach_segment(smap, out.origin, [&](uint_t i, vec1u ids) {
                if (ids.size() < params.min_area) {
                    smap[ids] = 0;
                    eids.push_back(i);
                }
            });

            inplace_remove(out.id, eids);
            inplace_remove(out.origin, eids);
            inplace_remove(out.px, eids);
            inplace_remove(out.py, eids);
            inplace_remove(out.area, eids);
        }

        return smap;
    }

    template <typename T, typename enable = typename std::enable_if<!std::is_pointer<T>::value>::type>
    vec2u segment(vec<2,T> map) {
        segment_output sdo; segment_params sdp;
        return segment(std::move(map), sdo, sdp);
    }

    struct segment_deblend_params {
        // Threshold value below which pixels will not be segmented
        double detect_threshold = 2.5;
        // Threshold distance to match a pixel to an existing segment (in pixels)
        double deblend_threshold = 5.0;
        // Minimum number of pixels per segment
        uint_t min_area = 0u;
        // First ID used to place segments on the map
        uint_t first_id = 1u;
        // Flux limit to compute the centroid (relative to peak; 0: all pixels in segmentation)
        double centroid_flim = 0.0;
    };

    struct segment_deblend_output {
        // ID of the segment as given on the segmentation map
        vec1u id;
        // Size of a segment in number of pixels
        vec1u area;
        // Barycenter
        vec1d bx, by;
        // Peak
        vec1i px, py;
        // Flat index of the peak value of a segment
        vec1u origin;
        // Total "flux" inside a segment
        vec1f flux;
    };

    inline vec2u segment_deblend(const vec2d& img, segment_deblend_output& out,
        const segment_deblend_params& params = segment_deblend_params()) {

        vec2u seg(img.dims);

        vif_check(params.first_id > 0, "first ID must be > 0");

        vec2u visited(img.dims);
        vec1u ox, oy;

        vec1u sid = sort(img);
        uint_t ipos = sid.size()-1;
        while (!is_finite(img.safe[sid.safe[ipos]])) {
            if (ipos == 0) break;

            --ipos;
        }

        uint_t mid = sid[ipos];

        auto explore_neighbors = [&](uint_t ty, uint_t tx) {
            auto check_add = [&](uint_t tty, uint_t ttx) {
                if (visited.safe(tty,ttx) != ipos) {
                    visited.safe(tty,ttx) = ipos;
                    ox.push_back(ttx);
                    oy.push_back(tty);
                }
            };

            if (ty != 0)             check_add(ty-1,tx);
            if (ty != img.dims[0]-1) check_add(ty+1,tx);
            if (tx != 0)             check_add(ty,tx-1);
            if (tx != img.dims[1]-1) check_add(ty,tx+1);
        };

        out.area.clear();
        out.origin.clear();
        out.id.data.reserve(0.1*sqrt(img.size()));
        out.area.data.reserve(0.1*sqrt(img.size()));
        out.origin.data.reserve(0.1*sqrt(img.size()));
        out.bx.data.reserve(0.1*sqrt(img.size()));
        out.by.data.reserve(0.1*sqrt(img.size()));
        out.px.data.reserve(0.1*sqrt(img.size()));
        out.py.data.reserve(0.1*sqrt(img.size()));

        vec1f flx_centroid;
        flx_centroid.data.reserve(0.1*sqrt(img.size()));

        uint_t id = params.first_id;

        while (img[mid] >= params.detect_threshold && ipos > 0) {
            // Check neighboring pixels within the deblend threshold if one of them
            // already has been segmented, and if so give this pixel to that segment
            vec1u tmid = mult_ids(img, mid);
            ox.clear(); oy.clear();
            oy.push_back(tmid[0]);
            ox.push_back(tmid[1]);
            visited[mid] = ipos;

            uint_t neib_seg = 0;
            double closest = dinf;
            uint_t iter = 0;
            while (!ox.empty() && neib_seg == 0 && iter < params.deblend_threshold) {
                vec1u tox = std::move(ox), toy = std::move(oy);
                ox.clear(); oy.clear();
                for (uint_t i : range(tox)) {
                    uint_t tseg = seg.safe(toy.safe[i],tox.safe[i]);
                    if (tseg > 0) {
                        // Favour the closest peak
                        vec1d ttid = mult_ids(img, out.origin.safe[tseg-params.first_id]);
                        double d = sqr(toy.safe[i] - ttid[0]) + sqr(tox.safe[i] - ttid[1]);
                        if (d < closest) {
                            neib_seg = tseg;
                            closest = d;
                        }
                    }

                    if (neib_seg == 0) {
                        explore_neighbors(toy.safe[i], tox.safe[i]);
                    }
                }

                ++iter;
            }

            if (neib_seg > 0) {
                seg[mid] = neib_seg;
            } else {
                seg[mid] = id;
                out.id.push_back(id);
                out.origin.push_back(mid);
                out.py.push_back(tmid[0]);
                out.px.push_back(tmid[1]);

                ++id;
            }

            --ipos;
            mid = sid[ipos];
        }

        // Erase islands
        {
            visited[_] = 0;
            vec2u old_seg = seg;

            seg[_] = 0;

            out.area.resize(out.id.size());
            out.flux.resize(out.id.size());
            out.by.resize(out.id.size());
            out.bx.resize(out.id.size());
            flx_centroid.resize(out.id.size());

            uint_t curseg = 0;

            auto fetch_neighbors = [&](uint_t ty, uint_t tx) {
                auto check_add = [&](uint_t tty, uint_t ttx) {
                    if (visited.safe(tty,ttx) != curseg) {
                        visited.safe(tty,ttx) = curseg;
                        ox.push_back(ttx);
                        oy.push_back(tty);
                    }
                };

                if (ty != 0)             check_add(ty-1,tx);
                if (ty != img.dims[0]-1) check_add(ty+1,tx);
                if (tx != 0)             check_add(ty,tx-1);
                if (tx != img.dims[1]-1) check_add(ty,tx+1);
            };

            for (uint_t s : range(out.area)) {
                curseg = out.id[s];

                vec1u tmid = mult_ids(seg, out.origin[s]);
                oy.clear(); ox.clear();
                oy.push_back(tmid[0]);
                ox.push_back(tmid[1]);

                while (!ox.empty()) {
                    vec1u tox = std::move(ox), toy = std::move(oy);
                    ox.clear(); oy.clear();
                    for (uint_t i : range(tox)) {
                        if (old_seg.safe(toy.safe[i],tox.safe[i]) == curseg) {
                            mid = flat_id(seg, toy[i], tox[i]);
                            seg[mid] = out.id[s];
                            out.area[s] += 1;
                            out.flux[s] += img[mid];
                            if (img[mid] >= params.centroid_flim*img[out.origin[s]]) {
                                out.by[s] += toy[i]*img[mid];
                                out.bx[s] += tox[i]*img[mid];
                                flx_centroid[s] += img[mid];
                            }

                            fetch_neighbors(toy[i], tox[i]);
                        }
                    }
                }
            }
        }

        // Erase too small regions
        if (count(out.area < params.min_area) != 0) {
            visited[_] = 0;

            uint_t curseg = 0;

            auto fetch_neighbors = [&](uint_t ty, uint_t tx) {
                auto check_add = [&](uint_t tty, uint_t ttx) {
                    if (visited.safe(tty,ttx) != curseg) {
                        visited.safe(tty,ttx) = curseg;
                        ox.push_back(ttx);
                        oy.push_back(tty);
                    }
                };

                if (ty != 0)             check_add(ty-1,tx);
                if (ty != img.dims[0]-1) check_add(ty+1,tx);
                if (tx != 0)             check_add(ty,tx-1);
                if (tx != img.dims[1]-1) check_add(ty,tx+1);
            };

            vec1u eids;
            for (uint_t s : range(out.area)) {
                curseg = out.id[s];

                if (out.area[s] >= params.min_area) continue;

                eids.push_back(s);

                vec1u tmid = mult_ids(img, out.origin[s]);
                ox.clear(); oy.clear();
                oy.push_back(tmid[0]);
                ox.push_back(tmid[1]);

                vec1u sids;
                sids.reserve(out.area[s]);

                uint_t neib_seg = 0;
                double closest = dinf;
                while (!ox.empty()) {
                    vec1u tox = std::move(ox), toy = std::move(oy);
                    ox.clear(); oy.clear();
                    for (uint_t i : range(tox)) {
                        uint_t tseg = seg.safe(toy.safe[i],tox.safe[i]);
                        if (tseg > 0) {
                            if (tseg != curseg) {
                                if (out.area.safe[tseg-params.first_id] >= params.min_area) {
                                    // Favour the closest peak
                                    vec1d ttid = mult_ids(img, out.origin.safe[tseg-params.first_id]);
                                    double d = sqr(toy.safe[i] - ttid[0]) + sqr(tox.safe[i] - ttid[1]);
                                    if (d < closest) {
                                        neib_seg = tseg;
                                        closest = d;
                                    }
                                }
                            } else {
                                sids.push_back(flat_id(img, toy.safe[i], tox.safe[i]));
                                fetch_neighbors(toy.safe[i], tox.safe[i]);
                            }
                        }
                    }
                }

                if (sids.size() != out.area[s]) {
                    error("huhuhh!!! ", sids.size(), " vs ", out.area[s]);
                }

                seg[sids] = neib_seg;
                if (neib_seg != 0) {
                    out.area[neib_seg-params.first_id] += out.area[s];
                    out.flux[neib_seg-params.first_id] += out.flux[s];
                    for (uint_t i : sids) {
                        if (img[i] >= params.centroid_flim*img[out.origin[neib_seg-params.first_id]]) {
                            tmid = mult_ids(img, i);
                            out.by[neib_seg-params.first_id] += tmid[0]*img[i];
                            out.bx[neib_seg-params.first_id] += tmid[1]*img[i];
                            flx_centroid[neib_seg-params.first_id] += img[i];
                        }
                    }
                }
            }

            inplace_remove(out.id, eids);
            inplace_remove(out.area, eids);
            inplace_remove(out.origin, eids);
            inplace_remove(flx_centroid, eids);
            inplace_remove(out.flux, eids);
            inplace_remove(out.by, eids);
            inplace_remove(out.bx, eids);
            inplace_remove(out.py, eids);
            inplace_remove(out.px, eids);
        }

        // Flux weighted barycenter
        out.by /= flx_centroid;
        out.bx /= flx_centroid;

        return seg;
    }

    inline void segment_distance(vec2u& map, vec2d& dmap, vec2u& imap) {
        dmap = replicate(dinf, map.dims);

        struct obj_state {
            uint_t id;
            uint_t npix = 0;
            std::vector<uint_t> oy, ox;
        };

        std::vector<obj_state> states;

        // Initialize states: identify segments and their boundaries where growth is allowed
        std::vector<uint_t> toy, tox;
        imap = replicate(npos, map.dims);
        vec2u omap = map;

        for (uint_t y : range(map.dims[0]))
        for (uint_t x : range(map.dims[1])) {
            if (omap.safe(y,x) == 0) {
                continue;
            }

            // Found a guy
            states.push_back(obj_state());
            auto& state = states.back();
            state.id = map.safe(y,x);

            toy.clear(); tox.clear();

            auto process_point = [&toy,&tox,&state,&map,&omap,&dmap,&imap](uint_t ty, uint_t tx) {
                omap.safe(ty,tx) = 0; // set to zero to avoid coming back to it
                dmap.safe(ty,tx) = 0;
                imap.safe(ty,tx) = flat_id(map, ty, tx);
                ++state.npix;

                auto check_add = [&toy,&tox,&state,&map,&omap,&dmap,&imap,ty,tx](uint_t tty, uint_t ttx) {
                    if (omap.safe(tty,ttx) == state.id) {
                        toy.push_back(tty);
                        tox.push_back(ttx);
                    } else if (map.safe(tty,ttx) == 0) {
                        state.oy.push_back(tty);
                        state.ox.push_back(ttx);
                        map.safe(tty,ttx) = state.id;
                        dmap.safe(tty,ttx) = 1.0;
                        imap.safe(tty,ttx) = flat_id(map, ty, tx);
                    }
                };

                if (ty != 0)             check_add(ty-1,tx);
                if (ty != map.dims[0]-1) check_add(ty+1,tx);
                if (tx != 0)             check_add(ty,tx-1);
                if (tx != map.dims[1]-1) check_add(ty,tx+1);
            };

            process_point(y, x);

            while (!tox.empty()) {
                uint_t ty = toy.back(); toy.pop_back();
                uint_t tx = tox.back(); tox.pop_back();
                process_point(ty, tx);
            }
        }

        // Now grow each segment one pixel at a time, and only keep the nearest in case of overlap
        bool starved = false;
        while (!starved) {
            starved = true;
            for (uint_t i : range(states)) {
                auto& state = states[i];

                if (state.ox.empty()) continue;
                starved = false;

                auto process_point = [&state,&map,&dmap,&imap](uint_t ty, uint_t tx) {
                    ++state.npix;

                    auto check_add = [&state,&map,&dmap,&imap,tx,ty](uint_t tty, uint_t ttx) {
                        double ox = imap.safe(ty,tx) % imap.dims[1];
                        double oy = imap.safe(ty,tx) / imap.dims[1];
                        double nd = sqr(double(tty) - oy) + sqr(double(ttx) - ox);

                        if (dmap.safe(tty,ttx) > nd) {
                            map.safe(tty,ttx) = state.id;
                            imap.safe(tty,ttx) = imap.safe(ty,tx);
                            dmap.safe(tty,ttx) = nd;
                            state.oy.push_back(tty);
                            state.ox.push_back(ttx);
                        }
                    };

                    if (ty != 0)             check_add(ty-1,tx);
                    if (ty != map.dims[0]-1) check_add(ty+1,tx);
                    if (tx != 0)             check_add(ty,tx-1);
                    if (tx != map.dims[1]-1) check_add(ty,tx+1);
                };

                toy = state.oy;   tox = state.ox;
                state.oy.clear(); state.ox.clear();

                while (!tox.empty()) {
                    uint_t ty = toy.back(); toy.pop_back();
                    uint_t tx = tox.back(); tox.pop_back();
                    process_point(ty, tx);
                }
            }
        }

        // map now contains the expanded segmentation map
        // imap now contains the pixel ID of the nearest neighbor
        // dmap now contains the squared distance to the nearest segment
    }
}

namespace impl {
    namespace astro_impl {
        template<typename T>
        struct extract_default_value {
            static constexpr const T value = 0;
        };

        template<>
        struct extract_default_value<float> {
            static constexpr const float value = fnan;
        };

        template<>
        struct extract_default_value<double> {
            static constexpr const double value = dnan;
        };

        template<>
        struct extract_default_value<std::string> {
            static constexpr const char* value = "";
        };
    }
}

namespace astro {
#ifndef NO_WCSLIB
    struct cutout_extractor {
        struct image_t {
            explicit image_t(const std::string& filename) :
                img(filename), hdr(img.read_header()), w(hdr), dims(img.image_dims()) {}

            image_t(image_t&& i) noexcept :
                img(std::move(i.img)), hdr(std::move(i.hdr)), w(std::move(i.w)),
                dims(std::move(i.dims)) {}

            fits::input_image img;
            fits::header      hdr;
            astro::wcs        w;
            vec1u             dims;
        };

        vec<1,image_t> imgs;
        double aspix = dnan;
        std::unique_ptr<image_t> dist;

        cutout_extractor() = default;

        void setup_image(const std::string& filename) {
            vec1s files;
            if (ends_with(filename, ".fits")) {
                files.push_back(filename);
            } else {
                files = fits::read_sectfits(filename);
            }

            imgs.data.reserve(files.size());
            imgs.dims[0] = files.size();
            for (auto& f : files) {
                imgs.data.emplace_back(f);
                if (!is_finite(aspix)) {
                    get_pixel_size(imgs.data.back().w, aspix);
                }
            }
        }

        void setup_distortion(const std::string& filename) {
            dist = std::unique_ptr<image_t>(new image_t(filename));
        }

    private:
        template<typename Type>
        bool get_cutout_(vec<2,Type>& cut, fits::header& hdr, bool gethdr,
            double ra, double dec, double size, Type def) const {

            cut.clear();

            double ira = ra, idec = dec;
            if (dist) {
                double dx, dy;
                astro::ad2xy(dist->w, ra, dec, dx, dy);
                dx -= 1.0; dy -= 1.0;
                if (dx > -0.5 && dy > -0.5 && dx < dist->dims[1]-0.5 && dy < dist->dims[0]-0.5) {
                    uint_t ix = round(dx), iy = round(dy);
                    dist->img.reach_hdu(1);
                    ira  -= dist->img.read_pixel({iy, ix});
                    dist->img.reach_hdu(2);
                    idec -= dist->img.read_pixel({iy, ix});
                }
            }

            int_t hs = max(1, ceil(0.5*size/aspix));
            double px = hs+1.0, py = hs+1.0;

            for (uint_t i : range(imgs)) {
                auto& iimg = imgs[i].img;
                auto& w = imgs[i].w;
                auto& dims = imgs[i].dims;

                double dxc, dyc;
                astro::ad2xy(w, ira, idec, dxc, dyc);
                dxc -= 1.0; dyc -= 1.0;
                int_t xc = round(dxc);
                int_t yc = round(dyc);
                int_t iy0 = yc-hs, iy1 = yc+hs, ix0 = xc-hs, ix1 = xc+hs;

                if (ix1 < 0 || ix0 >= int_t(dims[1]) || iy1 < 0 || iy0 >= int_t(dims[0])) {
                    // Source not covered
                    continue;
                }

                vec<2,Type> data;
                if (ix0 < 0 || ix1 >= int_t(dims[1]) || iy0 < 0 || iy1 >= int_t(dims[0])) {
                    // Source partially covered
                    data = replicate(def, 2*hs+1, 2*hs+1);

                    uint_t tx0 = max(0, ix0);
                    uint_t ty0 = max(0, iy0);
                    uint_t tx1 = min(dims[1]-1, ix1);
                    uint_t ty1 = min(dims[0]-1, iy1);

                    vec<2,Type> subcut;
                    iimg.read_subset(subcut, ty0-_-ty1, tx0-_-tx1);

                    uint_t x0 = int_t(tx0)-ix0, x1 = int_t(tx1)-ix0, y0 = int_t(ty0)-iy0, y1 = int_t(ty1)-iy0;
                    data(y0-_-y1,x0-_-x1) = std::move(subcut);
                } else {
                    // Source fully covered
                    int_t y0 = iy0, y1 = iy1, x0 = ix0, x1 = ix1;
                    iimg.read_subset(data, y0-_-y1, x0-_-x1);
                }

                if (cut.empty()) {
                    // First image on which this source is found
                    cut = std::move(data);

                    // Initialize WCS
                    if (gethdr) {
                        hdr = astro::filter_wcs(imgs[i].hdr);
                    }

                    // Save precise center
                    px = hs+1+(dxc-xc);
                    py = hs+1+(dyc-yc);
                } else {
                    // This source was found on another image, combine the two
                    vec1u idb = where(!is_finite(cut));
                    cut[idb] = data[idb];
                }
            }

            if (gethdr) {
                if (hdr.empty() && gethdr) {
                    // Source not covered, just use the first WCS header
                    hdr = astro::filter_wcs(imgs[0].hdr);
                }

                // Set cutout WCS
                if (!fits::setkey(hdr, "NAXIS1", cut.dims[1]) ||
                    !fits::setkey(hdr, "NAXIS2", cut.dims[0]) ||
                    !fits::setkey(hdr, "CRPIX1", px) || !fits::setkey(hdr, "CRPIX2", py) ||
                    !fits::setkey(hdr, "CRVAL1", ra) || !fits::setkey(hdr, "CRVAL2", dec)) {
                    return false;
                }
            }

            if (cut.empty()) {
                // Source not covered, just use placeholder data
                cut = replicate(def, 2*hs+1, 2*hs+1);
                return false;
            }

            return true;
        }

    public:
        template<typename Type, typename TypeD = Type>
        bool get_cutout(vec<2,Type>& cut, double ra, double dec, double size,
            TypeD def = impl::astro_impl::extract_default_value<Type>::value) const {
            fits::header hdr;
            return get_cutout_(cut, hdr, false, ra, dec, size, def);
        }

        template<typename Type, typename TypeD = Type>
        bool get_cutout(vec<2,Type>& cut, fits::header& hdr, double ra, double dec, double size,
            TypeD def = impl::astro_impl::extract_default_value<Type>::value) const {
            return get_cutout_(cut, hdr, true, ra, dec, size, def);
        }
    };
#else
    struct cutout_extractor {
        template<typename Dummy>
        void setup_image(const std::string&) {
            static_assert(!std::is_same<Dummy,Dummy>::value, "WCS support is disabled, "
                "please enable the WCSLib library to use this function");
        }

        template<typename Dummy>
        void setup_distortion(const std::string&) {
            static_assert(!std::is_same<Dummy,Dummy>::value, "WCS support is disabled, "
                "please enable the WCSLib library to use this function");
        }

        template<typename Type, typename TypeD = Type>
        bool get_cutout(vec<2,Type>&, double, double, double,
            TypeD def = impl::astro_impl::extract_default_value<Type>::value) const {
            static_assert(!std::is_same<Type,Type>::value, "WCS support is disabled, "
                "please enable the WCSLib library to use this function");

            return false;
        }

        template<typename Type, typename TypeD = Type>
        bool get_cutout(vec<2,Type>&, fits::header&, double, double, double,
            TypeD def = impl::astro_impl::extract_default_value<Type>::value) const {
            static_assert(!std::is_same<Type,Type>::value, "WCS support is disabled, "
                "please enable the WCSLib library to use this function");

            return false;
        }
    };
#endif
}

namespace impl {
    namespace astro_impl {
        // Convenience functions
        inline void regrid_drizzle_getbounds(const vec1d& v, uint_t& mi, int_t& ma, int_t n) {
            auto p = minmax(v);
            int_t tmp = floor(p.second+0.5);
            ma = tmp > n ? n : tmp;
            tmp = floor(p.first+0.5);
            mi = tmp > 0 ? tmp : 0;
        };

        // Function to find if a point lies on the right side of a polygon's edge
        // orient: orientation of the polygon
        // cx1, cy1, cx2, cy2: X and Y coordinates of the two nodes of the edge
        // x, y: X and Y coordinates of the point to test
        // return: true if the point is on the right side, false otherwise
        inline bool regrid_drizzle_in_poly_edge(int_t orient,
            double cx1, double cy1, double cx2, double cy2, double x, double y) {
            double cross = (cx2 - cx1)*(y - cy1) - (cy2 - cy1)*(x - cx1);
            return cross*orient < 0;
        };

        // Function to find the position and existence of the intersection of two lines
        // l1x1, l1y1, l1x2, l1y2: X and Y coordinates of two points defining the first line
        // l2x1, l2y1, l2x2, l2y2: X and Y coordinates of two points defining the second line
        // x, y: output X and Y coordinates of the intersection point (if any)
        // return: true if intersection exists, false otherwise
        inline bool regrid_drizzle_segment_intersect(
            double l1x1, double l1y1, double l1x2, double l1y2,
            double l2x1, double l2y1, double l2x2, double l2y2, double& x, double& y) {

            // Find the intersection point
            // http://stackoverflow.com/a/1968345/1565581
            double s1x = l1x2 - l1x1;
            double s1y = l1y2 - l1y1;
            double s2x = l2x2 - l2x1;
            double s2y = l2y2 - l2y1;

            double det = s1x*s2y - s1y*s2x;
            if (abs(det) < 5*std::numeric_limits<double>::epsilon()) {
                // 'det' is zero: the lines are parallel, no intersection
                return false;
            }

            double s12x = l1x1 - l2x1;
            double s12y = l1y1 - l2y1;
            double t = (s2x*s12y - s2y*s12x)/det;
            x = l1x1 + t*s1x;
            y = l1y1 + t*s1y;

            return true;
        };

        inline double polyon_area(const vec1d& x, const vec1d& y) {
            // vif_check(x.size() == y.size(), "incompatible dimensions between X and Y arrays (",
            //    x.size(), " vs. ", y.size(), ")");

            // Note: the above check is not performed here for performance reasons
            // and because this is an internal function. If you want to use it outside,
            // please uncomment the check

            double area = 0.0;

            if (x.size() < 3) return area;

            uint_t i1 = 0;
            for (uint_t i2 : range(1, x.size()-1)) {
                uint_t i3 = i2+1;

                area += 0.5*abs(
                    x.safe[i1]*(y.safe[i2]-y.safe[i3]) +
                    x.safe[i2]*(y.safe[i3]-y.safe[i1]) +
                    x.safe[i3]*(y.safe[i1]-y.safe[i2])
                );
            }

            return area;
        }

        template<typename Type>
        void add_flux(const vec<2,Type>& imgs, uint_t ix, uint_t iy, double flux_frac,
            vec<2,double>& res, uint_t idx, uint_t idy) {

            res.safe(idy,idx) += imgs.safe(iy,ix)*flux_frac;
        }

        template<typename Type>
        void add_flux(const vec<3,Type>& imgs, uint_t ix, uint_t iy, double flux_frac,
            vec<3,double>& res, uint_t idx, uint_t idy) {

            res.safe(_,idy,idx) += imgs.safe(_,iy,ix)*flux_frac;
        }

        template<std::size_t D, typename Type>
        inline void regrid_drizzle(const vec<D,Type>& imgs, uint_t ix, uint_t iy, double flux_frac,
            const vec1d& xps, const vec1d& yps, double pixel_area, vec<D,double>& res, vec2d& wei) {

            // Get bounds of this projection
            uint_t ymin, xmin;
            int_t tymax, txmax;
            regrid_drizzle_getbounds(yps, ymin, tymax, wei.dims[0]-1);
            regrid_drizzle_getbounds(xps, xmin, txmax, wei.dims[1]-1);

            if (tymax < 0 || txmax < 0 || xmin >= wei.dims[1] || ymin >= wei.dims[0]) {
                // Original pixel is outside of destination image
                return;
            }

            uint_t ymax = tymax, xmax = txmax;

            if (xmax == xmin && ymax == ymin) {
                // Original pixel is fully contained within destination pixel
                add_flux(imgs, ix, iy, flux_frac, res, xmin, ymin);
                wei.safe(ymin,xmin) += pixel_area;
                return;
            }

            // Get polygon orientation
            int_t orient = ((xps.safe[1] - xps.safe[0])*(yps.safe[2] - yps.safe[1]) -
                (yps.safe[1] - yps.safe[0])*(xps.safe[2] - xps.safe[1]) > 0) ? 1 : -1;

            // Sum flux from original pixel that fall inside the projection
            for (uint_t ipy = ymin; ipy <= ymax; ++ipy)
            for (uint_t ipx = xmin; ipx <= xmax; ++ipx) {
                // Construct the intersection polygon of the original pixel and the projection
                // https://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
                vec1d cpx = {ipx-0.5, ipx+0.5, ipx+0.5, ipx-0.5};
                vec1d cpy = {ipy-0.5, ipy-0.5, ipy+0.5, ipy+0.5};

                uint_t c2 = xps.size()-1;
                for (uint_t c1 : range(xps)) {
                    if (cpx.empty()) break;

                    vec1d icpx = cpx, icpy = cpy;
                    cpx.clear(); cpy.clear();

                    // Find out if all the current polygon's points are "inside" the
                    // projection's current edge
                    vec1b in(icpx.dims);
                    for (uint_t i : range(icpx)) {
                        in.safe[i] = regrid_drizzle_in_poly_edge(orient,
                            xps.safe[c1], yps.safe[c1], xps.safe[c2], yps.safe[c2],
                            icpx.safe[i], icpy.safe[i]
                        );
                    }

                    uint_t i2 = icpx.size()-1;
                    for (uint_t i1 : range(icpx)) {
                        if (in.safe[i2] != in.safe[i1]) {
                            // This edge [i2-i1] is intersected by the projection's
                            // current edge, find the intersection point and add it
                            // to the polygon
                            double tx, ty;
                            bool intersect = regrid_drizzle_segment_intersect(
                                xps.safe[c1], yps.safe[c1], xps.safe[c2], yps.safe[c2],
                                icpx.safe[i1], icpy.safe[i1], icpx.safe[i2], icpy.safe[i2],
                                tx, ty);

                            if (intersect) {
                                cpx.push_back(tx);
                                cpy.push_back(ty);
                            }
                        }

                        if (in.safe[i1]) {
                            // The point i1 is "inside" the projection's current edge,
                            // keep it for now
                            cpx.push_back(icpx.safe[i1]);
                            cpy.push_back(icpy.safe[i1]);
                        }

                        i2 = i1;
                    }

                    c2 = c1;
                }

                // No intersection, just discard that pixel
                if (cpx.size() < 3) continue;

                // Compute the area of this intersection
                double w = polyon_area(cpx, cpy);
                add_flux(imgs, ix, iy, flux_frac*(w/pixel_area), res, ipx, ipy);
                wei.safe(ipy,ipx) += w;
            }
        }


        template<typename T, typename U>
        bool regrid_nearest(const vec<2,T>& imgs, double xps, double yps, U& flx) {
            int_t mx = round(xps), my = round(yps);
            if (mx > 0 && uint_t(mx) < imgs.dims[1] && my > 0 && uint_t(my) < imgs.dims[0]) {
                flx = imgs.safe(uint_t(my), uint_t(mx));
                return true;
            } else {
                return false;
            }
        }

        template<typename T, typename U>
        bool regrid_linear(const vec<2,T>& imgs, double xps, double yps, U& flx) {
            flx = bilinear_strict(imgs, yps, xps, dnan);
            return is_finite(flx);
        }

        template<typename T, typename U>
        bool regrid_cubic(const vec<2,T>& imgs, double xps, double yps, U& flx) {
            flx = bicubic_strict(imgs, yps, xps, dnan);
            return is_finite(flx);
        }

        template<typename Type>
        void copy_subset(const vec<3,Type>& imgs, vec<3,double>& res, const vec1u& ids, const vec1u& idd) {
            uint_t nspix = imgs.dims[1]*imgs.dims[2];
            uint_t ndpix = res.dims[1]*res.dims[2];
            for (uint_t i : range(imgs.dims[0])) {
                res.safe[idd+i*ndpix] = imgs.safe[ids+i*nspix];
            }
        }

        template<typename Type>
        void copy_subset(const vec<2,Type>& imgs, vec<2,double>& res, const vec1u& ids, const vec1u& idd) {
            res.safe[idd] = imgs.safe[ids];
        }

        inline void apply_weights(vec3d& res, const vec2d& wei) {
            for (uint_t i : range(res.dims[0])) {
                res.safe(i,_,_) /= wei;
            }
        }

        inline void apply_weights(vec2d& res, const vec2d& wei) {
            res /= wei;
        }
    }
}

#ifndef NO_CFITSIO

namespace astro {

    enum class interpolation_method {
        nearest, linear, cubic
    };

    struct regrid_interpolate_params {
        bool verbose = false;
        bool conserve_flux = false;
        interpolation_method method = interpolation_method::linear;
        bool linearize = false;    // approximation: assume WCS transform is linear
    };

    template<typename T = double>
    vec2d regrid_interpolate(const vec<2,T>& imgs, const astro::wcs& astros, const astro::wcs& astrod,
        regrid_interpolate_params opts = regrid_interpolate_params{}) {

        // Regridded image
        vec2d res;

#ifdef NO_WCSLIB
        static_assert(!std::is_same<T,T>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else

        res = replicate(dnan, astrod.dims[astrod.y_axis], astrod.dims[astrod.x_axis]);

        uint_t nx = astrod.dims[astros.x_axis];
        uint_t ny = astrod.dims[astros.y_axis];

        // Precompute the projection of the new pixel grid on the old
        // Note: for the horizontal pixel 'i' of line 'j' (x_i,y_j), the grid is:
        //   (pux,puy)[i]     (pux,puy)[i+1]    # y_(j+0.5)
        //   (plx,ply)[i]     (plx,ply)[i+1]    # y_(j-0.5)
        //   # x_(i-0.5)      # x_(i+0.5)
        // To avoid re-computing stuff, pux and puy are moved into plx and ply
        // on each 'y' iteration for reuse.

        auto s2d_wcs = [&](double sx, double sy, double& dx, double& dy) {
            double tra, tdec;
            astro::xy2ad(astrod, sx+1, sy+1, tra, tdec);
            astro::ad2xy(astros, tra, tdec, dx, dy);
            dx -= 1.0; dy -= 1.0;
        };

        // If "linearize" is set, build approximate transformation matrix from image center
        double csx = nx/2;
        double csy = ny/2;
        double cdx, cdy, cdx1, cdy1, cdx2, cdy2;
        double pixel_area;
        if (opts.linearize) {
            s2d_wcs(csx, csy,     cdx,  cdy);
            s2d_wcs(csx+1.0, csy, cdx1, cdy1);
            s2d_wcs(csx, csy+1.0, cdx2, cdy2);

            cdx1 -= cdx; cdy1 -= cdy; cdx2 -= cdx; cdy2 -= cdy;
            pixel_area = sqrt(sqr(cdx1) + sqr(cdy1))*sqrt(sqr(cdx2) + sqr(cdy2));

            if (abs(cdx-round(cdx)) < 1e-6 && abs(cdy-round(cdy)) < 1e-6 &&
                abs(cdx1 - 1.0) < 1e-6 && abs(cdy2 - 1.0) < 1e-6 &&
                abs(cdy1) < 1e-3 && abs(cdx2) < 1e-3) {

                // This is a simple integer translation, let's optimize this
                int_t dx = round(cdx);
                int_t dy = round(cdy);
                int_t wx1 = nx/2;
                int_t wx2 = nx-1 - wx1;
                int_t wy1 = ny/2;
                int_t wy2 = ny-1 - wy1;

                vec1u ids, idd;
                subregion(res, {dy-wy1, dx-wx1, dy+wy2, dx+wx2}, idd, ids);

                res[_] = dnan;
                impl::astro_impl::copy_subset(imgs, res, ids, idd);

                return res;
            }
        }

        auto s2d = [&](double sx, double sy, double& dx, double& dy) {
            if (opts.linearize) {
                dx = (sx-csx)*cdx1 + (sy-csy)*cdx2 + cdx;
                dy = (sx-csx)*cdy1 + (sy-csy)*cdy2 + cdy;
            } else {
                s2d_wcs(sx, sy, dx, dy);
            }
        };

        auto pg = progress_start(res.size());
        vec1d plx(nx+1);
        vec1d ply(nx+1);
        for (uint_t ix : range(nx+1)) {
            s2d(ix-0.5, -0.5, plx.safe[ix], ply.safe[ix]);
        }

        vec1d pux(nx+1);
        vec1d puy(nx+1);
        for (uint_t iy : range(ny)) {
            for (uint_t ix : range(nx+1)) {
                s2d(ix-0.5, iy+0.5, pux.safe[ix], puy.safe[ix]);
            }

            for (uint_t ix : range(nx)) {
                // Find projection of each pixel of the new grid on the original image
                // NB: assumes the astrometry is such that this projection is
                // reasonably approximated by a 4-edge polygon (i.e.: varying pixel scales,
                // pixel offsets and rotations are fine, but weird things may happen close
                // to the poles of the projection where things become non-linear)

                double xps = 0.25*(plx.safe[ix] + plx.safe[ix+1] + pux.safe[ix+1] + pux.safe[ix]);
                double yps = 0.25*(ply.safe[ix] + ply.safe[ix+1] + puy.safe[ix+1] + puy.safe[ix]);

                double flx = 0.0;
                bool covered = false;

                switch (opts.method) {
                case interpolation_method::nearest:
                    covered = impl::astro_impl::regrid_nearest(imgs, xps, yps, flx);
                    break;
                case interpolation_method::linear:
                    covered = impl::astro_impl::regrid_linear(imgs, xps, yps, flx);
                    break;
                case interpolation_method::cubic:
                    covered = impl::astro_impl::regrid_cubic(imgs, xps, yps, flx);
                    break;
                }

                if (covered) {
                    if (opts.conserve_flux) {
                        if (!opts.linearize) {
                            pixel_area = impl::astro_impl::polyon_area(
                                {plx.safe[ix], plx.safe[ix+1], pux.safe[ix+1], pux.safe[ix]},
                                {ply.safe[ix], ply.safe[ix+1], puy.safe[ix+1], puy.safe[ix]}
                            );
                        }

                        flx *= pixel_area;
                    }

                    res.safe(iy,ix) = flx;
                }

                if (opts.verbose) progress(pg, 31);
            }

            std::swap(plx, pux);
            std::swap(ply, puy);
        }
#endif

        return res;
    }

    struct regrid_drizzle_params {
        bool verbose = false;      // print progress in terminal
        double pixfrac = 1.0;      // fraction of pixel to drizzle (1 = simple projection)
        bool linearize = false;    // approximation: assume WCS transform is linear
        bool dest_pixfrac = false; // approximation: assume pixfrac is linear on destination grid
    };

    template<std::size_t D, typename T = double>
    vec<D,double> regrid_drizzle(const vec<D,T>& imgs, const astro::wcs& astros, const astro::wcs& astrod,
        vec2d& wei, regrid_drizzle_params opts = regrid_drizzle_params{}) {

        // Regridded image
        vec<D,double> res;

#ifdef NO_WCSLIB
        static_assert(!std::is_same<T,T>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else

        vif_check(D == astros.dims.size(), "mismatch between WCS and image number of dimensions (",
            astros.dims.size(), " vs. ", D, ")");

        for (uint_t i : range(D)) {
            vif_check(imgs.dims[i] == astros.dims[i], "mismatch between WCS and image dimension (",
                astros.dims, " vs. ", imgs.dims, ")");
        }

        auto resd = imgs.dims;
        resd[D-2] = astrod.dims[astrod.y_axis];
        resd[D-1] = astrod.dims[astrod.x_axis];
        res = replicate(0.0, resd);
        wei = replicate(0.0, astrod.dims[astrod.y_axis], astrod.dims[astrod.x_axis]);

        uint_t nx = astros.dims[astros.x_axis];
        uint_t ny = astros.dims[astros.y_axis];

        // Precompute the projection of the old pixel grid on the new
        // In case pixfrac=1:
        // For the horizontal pixel 'i' of line 'j' (x_i,y_j), the grid is:
        //   (pux,puy)[i]     (pux,puy)[i+1]    # y_(j+0.5)
        //   (plx,ply)[i]     (plx,ply)[i+1]    # y_(j-0.5)
        //   # x_(i-0.5)      # x_(i+0.5)
        // To avoid re-computing stuff, pux and puy are swapped into plx and ply
        // on each 'y' iteration for reuse. Total of (Nx+1)*(Ny+1) WCS transforms.
        //
        // In case pixfrac!=1:
        // If dest_pixfrac=true, then the projection's coordinates are
        // computed by linearly interpolating the above.
        // If dest_pixfrac=false:
        // For the horizontal pixel 'i' of line 'j' (x_i,y_j), the grid is:
        //   (plux,pluy)[i]   (prux,pruy)[i]    # y_(j+0.5*pf)
        //   (pllx,plly)[i]   (prlx,prly)[i]    # y_(j-0.5*pf)
        //   # x_(i-0.5*pf)   # x_(i+0.5*pf)
        // Total of 4*Nx*Ny WCS transforms.

        auto s2d_wcs = [&](double sx, double sy, double& dx, double& dy) {
            double tra, tdec;
            astro::xy2ad(astros, sx+1, sy+1, tra, tdec);
            astro::ad2xy(astrod, tra, tdec, dx, dy);
            dx -= 1.0; dy -= 1.0;
        };

        // If "linearize" is set, build approximate transformation matrix from image center
        double csx = nx/2;
        double csy = ny/2;
        double cdx, cdy, cdx1, cdy1, cdx2, cdy2;
        double pixel_area;
        if (opts.linearize) {
            s2d_wcs(csx, csy,     cdx,  cdy);
            s2d_wcs(csx+1.0, csy, cdx1, cdy1);
            s2d_wcs(csx, csy+1.0, cdx2, cdy2);

            cdx1 -= cdx; cdy1 -= cdy; cdx2 -= cdx; cdy2 -= cdy;
            pixel_area = sqrt(sqr(cdx1) + sqr(cdy1))*sqrt(sqr(cdx2) + sqr(cdy2));

            if (abs(cdx-round(cdx)) < 1e-6 && abs(cdy-round(cdy)) < 1e-6 &&
                abs(cdx1 - 1.0) < 1e-6 && abs(cdy2 - 1.0) < 1e-6 &&
                abs(cdy1) < 1e-3 && abs(cdx2) < 1e-3) {

                // This is a simple integer translation, let's optimize this
                int_t dx = round(cdx);
                int_t dy = round(cdy);
                int_t wx1 = nx/2;
                int_t wx2 = nx-1 - wx1;
                int_t wy1 = ny/2;
                int_t wy2 = ny-1 - wy1;

                vec1u ids, idd;
                subregion(wei, {dy-wy1, dx-wx1, dy+wy2, dx+wx2}, idd, ids);

                res[_] = dnan;
                impl::astro_impl::copy_subset(imgs, res, ids, idd);
                wei[idd] = 1.0;

                return res;
            }
        }

        auto s2d = [&](double sx, double sy, double& dx, double& dy) {
            if (opts.linearize) {
                dx = (sx-csx)*cdx1 + (sy-csy)*cdx2 + cdx;
                dy = (sx-csx)*cdy1 + (sy-csy)*cdy2 + cdy;
            } else {
                s2d_wcs(sx, sy, dx, dy);
            }
        };

        auto pg = progress_start(nx*ny);
        vec1d plx, ply;
        vec1d pux, puy;
        vec1d pllx, prlx, plux, prux;
        vec1d plly, prly, pluy, pruy;
        if (opts.pixfrac == 1 || opts.dest_pixfrac) {
            plx.resize(nx+1);
            ply.resize(nx+1);
            for (uint_t ix : range(nx+1)) {
                s2d(ix-0.5, -0.5, plx.safe[ix], ply.safe[ix]);
            }

            pux.resize(nx+1);
            puy.resize(nx+1);
        } else {
            pllx.resize(nx); prlx.resize(nx);
            plux.resize(nx); prux.resize(nx);
            plly.resize(nx); prly.resize(nx);
            pluy.resize(nx); pruy.resize(nx);
        }

        for (uint_t iy : range(ny)) {
            if (opts.pixfrac == 1 || opts.dest_pixfrac) {
                for (uint_t ix : range(nx+1)) {
                    s2d(ix-0.5, iy+0.5, pux.safe[ix], puy.safe[ix]);
                }
            } else {
                const double dp = 0.5*opts.pixfrac;
                for (uint_t ix : range(nx)) {
                    s2d(ix-dp, iy-dp, pllx.safe[ix], plly.safe[ix]);
                    s2d(ix+dp, iy-dp, prlx.safe[ix], prly.safe[ix]);
                    s2d(ix-dp, iy+dp, plux.safe[ix], pluy.safe[ix]);
                    s2d(ix+dp, iy+dp, prux.safe[ix], pruy.safe[ix]);
                }
            }

            for (uint_t ix : range(nx)) {
                // Find projection of each pixel of the original image on the new grid
                // NB: assumes the astrometry is such that this projection is
                // reasonably approximated by a 4-edge polygon (i.e.: varying pixel scales,
                // pixel offsets and rotations are fine, but weird things may happen close
                // to the poles of the projection where things become non-linear)

                vec1d xps, yps;
                if (opts.pixfrac == 1 || opts.dest_pixfrac) {
                    xps = {plx.safe[ix], plx.safe[ix+1], pux.safe[ix+1], pux.safe[ix]};
                    yps = {ply.safe[ix], ply.safe[ix+1], puy.safe[ix+1], puy.safe[ix]};

                    if (opts.dest_pixfrac && opts.pixfrac != 1) {
                        double mx = mean(xps);
                        xps = mx + opts.pixfrac*(xps - mx);
                        double my = mean(yps);
                        yps = my + opts.pixfrac*(yps - my);
                    }
                } else {
                    xps = {pllx.safe[ix], prlx.safe[ix], prux.safe[ix], plux.safe[ix]};
                    yps = {plly.safe[ix], prly.safe[ix], pruy.safe[ix], pluy.safe[ix]};
                }

                if (!opts.linearize) {
                    pixel_area = impl::astro_impl::polyon_area(xps, yps);
                }

                double flx_frac = sqr(opts.pixfrac);
                impl::astro_impl::regrid_drizzle(imgs, ix, iy, flx_frac, xps, yps, pixel_area, res, wei);

                if (opts.verbose) progress(pg, 31);
            }

            if (opts.pixfrac == 1 || opts.dest_pixfrac) {
                std::swap(plx, pux);
                std::swap(ply, puy);
            }
        }
#endif

        impl::astro_impl::apply_weights(res, wei);

        return res;
    }

    template<std::size_t D, typename T = double>
    vec<D,double> regrid_drizzle(const vec<D,T>& imgs, const astro::wcs& astros, const astro::wcs& astrod,
        regrid_drizzle_params opts = regrid_drizzle_params{}) {
        vec2d wei;
        return regrid_drizzle(imgs, astros, astrod, wei, opts);
    }
}

#endif

}

#endif
