#ifndef IMAGE_HPP
#define IMAGE_HPP

#include "vec.hpp"
#include "math.hpp"

template<typename Type>
typename vec_t<2,Type>::effective_type enlarge(const vec_t<2,Type>& v, int_t pix,
    const typename vec_t<2,Type>::rtype& def = 0.0) {

    if (pix >= 0) {
        uint_t upix = pix;

        typename vec_t<2,Type>::effective_type r = dblarr(v.dims[0]+2*upix, v.dims[1]+2*upix) + def;

        for (uint_t y : range(v.dims[0]))
        for (uint_t x : range(v.dims[1])) {
            r(x+upix,y+upix) = v(x,y);
        }

        return r;
    } else {
        uint_t upix = -pix;
        phypp_check(2*upix < v.dims[0] && 2*upix < v.dims[1],
            "cannot shrink image to negative size (", upix, " vs [", v.dims, "]");

        typename vec_t<2,Type>::effective_type r = dblarr(v.dims[0]-2*upix, v.dims[1]-2*upix);

        for (uint_t y : range(r.dims[0]))
        for (uint_t x : range(r.dims[1])) {
            r(x,y) = v(x+upix,y+upix);
        }

        return r;
    }
}

// Get a sub region 'reg' inside an image. reg = {x0, y0, x1, y1}.
template<typename TypeV, typename TypeR = int_t>
void subregion(const vec_t<2,TypeV>& v, const vec_t<1,TypeR>& reg, vec1u& rr, vec1u& rs) {
    assert(n_elements(reg) == 4);
    int_t nvx = v.dims[0], nvy = v.dims[1];
    int_t nx = reg[2]-reg[0]+1, ny = reg[3]-reg[1]+1;

    vec1i vreg = reg;
    vec1i sreg = {0,0,nx-1,ny-1};

    if (reg[0] >= nvx || reg[0] > reg[2] || reg[1] >= nvy || reg[1] > reg[3]) {
        rr = rs = uintarr(0);
        return;
    }

    if (reg[0] < 0) {
        vreg[0] = 0;
        sreg[0] += 0 - reg[0];
    }
    if (reg[2] >= nvx) {
        vreg[2] = nvx-1;
        sreg[2] -= reg[2] - (nvx-1);
    }

    if (reg[1] < 0) {
        vreg[1] = 0;
        sreg[1] += 0 - reg[1];
    }
    if (reg[3] >= nvy) {
        vreg[3] = nvy-1;
        sreg[3] -= reg[3] - (nvy-1);
    }

    vec1u vx = rgen(vreg[0], vreg[2]);
    vec1u vy = rgen(vreg[1], vreg[3]);
    vec1u sx = rgen(sreg[0], sreg[2]);
    vec1u sy = rgen(sreg[1], sreg[3]);

    rr = flatten(uindgen(nvx,nvy)(vx,vy));
    rs = flatten(uindgen(nx,ny)(sx,sy));
}

template<typename TypeV, typename TypeR>
typename vec_t<2,TypeV>::effective_type subregion(const vec_t<2,TypeV>& v,
    const vec_t<1,TypeR>& reg, const typename vec_t<2,TypeV>::rtype& def = 0.0) {

    vec1u rr, rs;
    subregion(v, reg, rr, rs);

    int_t nx = reg(2)-reg(0)+1, ny = reg(3)-reg(1)+1;
    vec_t<2,rtype_t<TypeV>> sub = replicate(rtype_t<TypeV>(def), nx, ny);

    sub[rs] = v[rr];

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

        trs(x,y) = v(rx,ry)*(nx - tx)*(ny - ty) +
                   v(nx,ry)*(tx - rx)*(ny - ty) +
                   v(rx,ny)*(nx - tx)*(ty - ry) +
                   v(nx,ny)*(tx - rx)*(ty - ry);
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

    vec3d m = dblarr(dim[0], dim[1], 4);

    vec2d px = replicate(dindgen(dim[1]), dim[0]);
    vec2d py = transpose(replicate(dindgen(dim[0]), dim[1]));

    m(_,_,0) = pow(px-0.5 - center[1],2) + pow(py-0.5 - center[0],2) <= radius*radius;
    m(_,_,1) = pow(px+0.5 - center[1],2) + pow(py-0.5 - center[0],2) <= radius*radius;
    m(_,_,2) = pow(px+0.5 - center[1],2) + pow(py+0.5 - center[0],2) <= radius*radius;
    m(_,_,3) = pow(px-0.5 - center[1],2) + pow(py+0.5 - center[0],2) <= radius*radius;

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
        res[i] = total(mask*img)/total(mask);
    }

    return res;
}

template<typename F>
auto generate_img(const vec1u& dims, F&& expr) -> vec_t<2,decltype(expr(0,0))> {
    phypp_check(dims.size() == 2 || dims.size() == 1,
        "generate_img: must provide one or two numbers for image size");

    vec_t<2,decltype(expr(0,0))> img;
    if (dims.size() == 1) img.resize(dims[0], dims[0]);
    else img.resize(dims[0], dims[1]);

    for (uint_t x : range(img.dims[0]))
    for (uint_t y : range(img.dims[1])) {
        img(x,y) = expr(x,y);
    }

    return img;
}

#endif
