#ifndef IMAGE_HPP
#define IMAGE_HPP

#include "phypp/vec.hpp"
#include "phypp/math.hpp"

template<typename Type>
vec_t<2,rtype_t<Type>> enlarge(const vec_t<2,Type>& v, const std::array<uint_t,4> upix,
    const rtype_t<Type>& def = 0.0) {

    vec_t<2,rtype_t<Type>> r = replicate(def,
        v.dims[0]+upix[0]+upix[2], v.dims[1]+upix[1]+upix[3]);

    for (uint_t x : range(v.dims[0]))
    for (uint_t y : range(v.dims[1])) {
        r.safe(x+upix[0],y+upix[1]) = v.safe(x,y);
    }

    return r;
}

template<typename Type>
vec_t<2,rtype_t<Type>> enlarge(const vec_t<2,Type>& v, uint_t upix,
    const rtype_t<Type>& def = 0.0) {
    return enlarge(v, {{upix, upix, upix, upix}}, def);
}

template<typename Type>
vec_t<2,rtype_t<Type>> shrink(const vec_t<2,Type>& v, const std::array<uint_t,4> upix,
    const rtype_t<Type>& def = 0.0) {

    if (upix[0]+upix[2] >= v.dims[0] || upix[1]+upix[3] >= v.dims[1]) {
        return vec_t<2,rtype_t<Type>>();
    }

    vec_t<2,rtype_t<Type>> r(v.dims[0]-upix[0]-upix[2], v.dims[1]-upix[1]-upix[3]);

    for (uint_t x : range(r.dims[0]))
    for (uint_t y : range(r.dims[1])) {
        r.safe(x,y) = v.safe(x+upix[0],y+upix[1]);
    }

    return r;
}

template<typename Type>
vec_t<2,rtype_t<Type>> shrink(const vec_t<2,Type>& v, uint_t upix,
    const rtype_t<Type>& def = 0.0) {

    return shrink(v, {{upix, upix, upix, upix}}, def);
}

// Get a sub region 'reg' inside an image. reg = {x0, y0, x1, y1}.
template<typename TypeV, typename TypeR = int_t>
void subregion(const vec_t<2,TypeV>& v, const vec_t<1,TypeR>& reg, vec1u& rr, vec1u& rs) {
    assert(n_elements(reg) == 4);
    int_t nvx = v.dims[0], nvy = v.dims[1];
    int_t nx = reg.safe[2]-reg.safe[0]+1, ny = reg.safe[3]-reg.safe[1]+1;

    vec1i vreg = reg;
    vec1i sreg = {0,0,nx-1,ny-1};

    if (reg.safe[0] >= nvx || reg.safe[0] > reg.safe[2] ||
        reg.safe[1] >= nvy || reg.safe[1] > reg.safe[3]) {
        rr = rs = uintarr(0);
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

    vec1u vx = rgen(vreg.safe[0], vreg.safe[2]);
    vec1u vy = rgen(vreg.safe[1], vreg.safe[3]);
    vec1u sx = rgen(sreg.safe[0], sreg.safe[2]);
    vec1u sy = rgen(sreg.safe[1], sreg.safe[3]);

    rr = flatten(uindgen(nvx,nvy).safe(vx,vy));
    rs = flatten(uindgen(nx,ny).safe(sx,sy));
}

template<typename TypeV, typename TypeR>
typename vec_t<2,TypeV>::effective_type subregion(const vec_t<2,TypeV>& v,
    const vec_t<1,TypeR>& reg, const typename vec_t<2,TypeV>::rtype& def = 0.0) {

    vec1u rr, rs;
    subregion(v, reg, rr, rs);

    int_t nx = reg(2)-reg(0)+1, ny = reg(3)-reg(1)+1;
    vec_t<2,rtype_t<TypeV>> sub = replicate(rtype_t<TypeV>(def), nx, ny);

    sub.safe[rs] = v.safe[rr];

    return sub;
}

template<typename TypeV, typename TypeR>
typename vec_t<2,TypeV>::effective_type subregion(const vec_t<2,TypeV>& v,
    std::initializer_list<TypeR>&& reg, const typename vec_t<2,TypeV>::rtype& def = 0.0) {

    return subregion(v, vec_t<1,TypeR>(std::move(reg)), def);
}

template<typename TypeV, typename TypeD = double>
typename vec_t<2,TypeV>::effective_type translate(const vec_t<2,TypeV>& v, double dx, double dy,
    const typename vec_t<2,TypeV>::rtype& def = 0.0) {

    vec_t<2,rtype_t<TypeV>> trs = replicate(rtype_t<TypeV>(def), v.dims);
    for (uint_t x : range(v.dims[0]))
    for (uint_t y : range(v.dims[1])) {
        double tx = x - dx;
        double ty = y - dy;
        int_t rx = round(tx);
        int_t ry = round(ty);

        if (rx < 0 || rx > int_t(v.dims[0])-1 || ry < 0 || ry > int_t(v.dims[1])-1) continue;
        int_t nx, ny;
        if (tx - rx > 0) {
            nx = rx + 1;
            if (nx > int_t(v.dims[0]-1)) continue;
        } else {
            nx = rx;
            rx = nx - 1;
            if (rx < 0) continue;
        }
        if (ty - ry > 0) {
            ny = ry + 1;
            if (ny > int_t(v.dims[1]-1)) continue;
        } else {
            ny = ry;
            ry = ny - 1;
            if (ry < 0) continue;
        }

        trs.safe(x,y) = v.safe(rx,ry)*(nx - tx)*(ny - ty) +
                        v.safe(nx,ry)*(tx - rx)*(ny - ty) +
                        v.safe(rx,ny)*(nx - tx)*(ty - ry) +
                        v.safe(nx,ny)*(tx - rx)*(ty - ry);
    }

    return trs;
}

vec2d circular_mask(vec1u dim, const vec1d& center, double radius) {
    if (n_elements(dim) == 0) {
        dim = {uint_t(radius), uint_t(radius)};
    } else if (n_elements(dim) == 1) {
        dim = {dim[0], dim[0]};
    }

    assert(n_elements(center) == 2);

    vec3d m(dim[0], dim[1], 4);

    vec2d px = replicate(dindgen(dim[1]), dim[0]);
    vec2d py = transpose(replicate(dindgen(dim[0]), dim[1]));

    m.safe(_,_,0) = pow(px-0.5 - center.safe[1],2) + pow(py-0.5 - center.safe[0],2) <= radius*radius;
    m.safe(_,_,1) = pow(px+0.5 - center.safe[1],2) + pow(py-0.5 - center.safe[0],2) <= radius*radius;
    m.safe(_,_,2) = pow(px+0.5 - center.safe[1],2) + pow(py+0.5 - center.safe[0],2) <= radius*radius;
    m.safe(_,_,3) = pow(px-0.5 - center.safe[1],2) + pow(py+0.5 - center.safe[0],2) <= radius*radius;

    return mean(m,2);
}

template<typename Type>
vec_t<1, rtype_t<Type>> radial_profile(const vec_t<2,Type>& img, uint_t npix) {
    vec_t<1, rtype_t<Type>> res(npix);
    uint_t hsx = img.dims[0]/2;
    uint_t hsy = img.dims[1]/2;
    res[0] = img(hsx,hsy);
    for (uint_t i : range(1u, npix)) {
        vec2d mask = circular_mask({img.dims[0], img.dims[1]}, {double(hsx), double(hsy)}, i)*
            (1.0 - circular_mask({img.dims[0], img.dims[1]}, {double(hsx), double(hsy)}, i-1));
        res.safe[i] = total(mask*img)/total(mask);
    }

    return res;
}

template<typename F>
auto generate_img(const std::array<uint_t,2>& dims, F&& expr) -> vec_t<2,decltype(expr(0,0))> {
    vec_t<2,decltype(expr(0,0))> img(dims);
    for (uint_t x : range(img.dims[0]))
    for (uint_t y : range(img.dims[1])) {
        img.safe(x,y) = expr(x,y);
    }

    return img;
}

vec2d gaussian_profile(const std::array<uint_t,2>& dims, double sigma) {
    return generate_img(dims, [&](int_t x, int_t y) {
        return exp(-(sqr(x - int_t(dims[0]/2)) + sqr(y - int_t(dims[1]/2)))/(2.0*sqr(sigma)));
    });
}

// Perform the convolution of two 2D arrays, assuming the second one is the kernel.
// Naive loop implementation: this is slow but simple and reliable.
template<typename TypeY1, typename TypeY2>
auto convolve2d_naive(const vec_t<2,TypeY1>& map, const vec_t<2,TypeY2>& kernel) ->
    vec_t<2,decltype(map[0]*kernel[0])> {
    phypp_check(kernel.dims[0]%2 == 1, "kernel must have odd dimensions (", kernel.dims, ")");
    phypp_check(kernel.dims[1]%2 == 1, "kernel must have odd dimensions (", kernel.dims, ")");

    uint_t hxsize = kernel.dims[0]/2;
    uint_t hysize = kernel.dims[1]/2;

    vec_t<2,decltype(map[0]*kernel[0])> r(map.dims);
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
auto convolve2d(const vec_t<2,TypeY1>& map, const vec_t<2,TypeY2>& kernel) ->
    vec_t<2,decltype(map[0]*kernel[0])> {
#ifdef NO_FFTW
    return convolve2d_naive(map, kernel);
#else
    phypp_check(kernel.dims[0]%2 == 1, "kernel must have odd dimensions (", kernel.dims, ")");
    phypp_check(kernel.dims[1]%2 == 1, "kernel must have odd dimensions (", kernel.dims, ")");

    uint_t hsx = kernel.dims[0]/2, hsy = kernel.dims[1]/2;

    // Pad image to prevent issues with cyclic borders
    vec2d tmap = enlarge(map, {{hsx, hsy, hsx, hsy}});

    // Resize kernel to map size, with kernel center at (0,0)
    vec2d tkernel(tmap.dims);

    vec1u px1 = hsx + uindgen(hsx);
    vec1u py1 = hsy + uindgen(hsy);
    vec1u px2 = hsx - 1 - uindgen(hsx);
    vec1u py2 = hsy - 1 - uindgen(hsy);

    vec1u ix1 = uindgen(hsx);
    vec1u iy1 = uindgen(hsy);
    vec1u ix2 = tmap.dims[0] - 1 - uindgen(hsx);
    vec1u iy2 = tmap.dims[1] - 1 - uindgen(hsy);

    tkernel(ix1,iy1) = kernel(px1, py1);
    tkernel(ix2,iy1) = kernel(px2, py1);
    tkernel(ix1,iy2) = kernel(px1, py2);
    tkernel(ix2,iy2) = kernel(px2, py2);

    // Perform the convolution in Fourrier space
    auto cimg = fft(tmap)*fft(tkernel);

    // Go back to real space and shrink map back to original dimensions
    return shrink(ifft(cimg), {{hsx, hsy, hsx, hsy}})/cimg.size();
#endif
}

template<typename T, typename F>
auto boxcar(const vec_t<2,T>& img, uint_t hsize, F&& func) ->
    vec_t<2,decltype(func(flatten(img)))> {

    vec_t<2,decltype(func(flatten(img)))> res(img.dims);
    for (uint_t x = 0; x < img.dims[0]; ++x)
    for (uint_t y = 0; y < img.dims[1]; ++y) {
        uint_t x0 = x >= hsize ? x - hsize : 0;
        uint_t x1 = x < img.dims[0]-hsize ? x + hsize : img.dims[0]-1;
        uint_t y0 = y >= hsize ? y - hsize : 0;
        uint_t y1 = y < img.dims[1]-hsize ? y + hsize : img.dims[1]-1;

        vec_t<1,rtype_t<T>> tmp;
        tmp.reserve((x1-x0)*(y1-y0));
        for (uint_t tx = x0; tx <= x1; ++tx)
        for (uint_t ty = y0; ty <= y1; ++ty) {
            tmp.push_back(img.safe(tx,ty));
        }

        res.safe(x,y) = func(tmp);
    }

    return res;
}

#endif
